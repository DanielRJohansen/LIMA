#include "Lipids.h"
#include "Format.h"
#include "Simulation.cuh"
#include "Statistics.h"
#include "MoleculeUtils.h"
#include "MoleculeGraph.h"

#include "Display.h"	// _MakeLipid shows the lipid

Lipids::Select::Select(const std::string& lipidname, const fs::path& workDir, double percentage) :
	userSupplied(fs::exists(workDir / (lipidname + ".gro")) && fs::exists(workDir / (lipidname + ".itp"))),
	lipidname(lipidname),
	percentage(percentage)
{
	const fs::path defaultLipidsDir = FileUtils::GetLimaDir() / ("resources/Slipids");

	if (!userSupplied && !fs::exists(defaultLipidsDir / (lipidname + ".itp"))) {
		throw std::runtime_error(Lima::Format("Failed to find lipid: {}, looked here: \n\t{}\nAnd here:\n\t{}",
			lipidname, workDir.string(), defaultLipidsDir.string()));
	}

	const fs::path& lipidDir = userSupplied ? workDir : defaultLipidsDir;
	grofile = std::make_unique<GroFile>(lipidDir / (lipidname + ".gro"));
	topfile = std::make_unique<TopologyFile>(lipidDir / (lipidname + ".itp"));

	if (userSupplied) 
		OrganizeLipidIntoCompoundsizedSections(*grofile, topfile->GetMoleculeType());
}

std::string Lipids::NameSelection(const Selection& selection) {
	std::ostringstream oss;
	for (const auto& sel : selection) {
		oss << sel.lipidname << "(" << Lima::Format("{:.1f}", sel.percentage) << "%)_";
	}
	return oss.str();
}

Lipids::Selection Lipids::Composition::Load() const {
	Selection selection;
	for (const Lipid& lipid : lipids)
		selection.emplace_back(lipid.name, workDir, lipid.percentage);
	return selection;
}

std::string Lipids::Composition::Name() const {
	std::ostringstream oss;
	for (const Lipid& lipid : lipids)
		oss << lipid.name << "(" << Lima::Format("{:.1f}", lipid.percentage) << "%)_";
	return oss.str();
}



Float3 LipidheadCenterOfMassEstimate(const GroFile& grofile) {
	Float3 sum{ 0 };
	int count = 0;
	for (const auto& atom : grofile.atoms) {
		if (atom.atomName[0] != 'C' && atom.atomName[0] != 'H') {
			sum += atom.position;
			count++;
		}
	}
	return sum / static_cast<float>(count);
}

Float3 LipidbodyGeometricCenter(const GroFile& grofile) {
	Float3 boundingboxMax{ -FLT_MAX }, boundingboxMin{ FLT_MAX };

	for (const auto& atom : grofile.atoms) {
		if (atom.atomName[0] == 'C') {
			boundingboxMax = Float3::ElementwiseMax(boundingboxMax, atom.position);
			boundingboxMin = Float3::ElementwiseMin(boundingboxMin, atom.position);
		}
	}

	return (boundingboxMax + boundingboxMin) / 2.f;
}


void Lipids::OrientLipidhead(GroFile& grofile, Float3 desiredOrientation) {
	// These are NOT precise, just naive heuristics
	const Float3 headCenter = LipidheadCenterOfMassEstimate(grofile);
	const Float3 lipidbodyCenter = LipidbodyGeometricCenter(grofile);

	if ((headCenter- lipidbodyCenter).len() < 0.1f) {
		// The estimate of the head relative to the body is too weak to orient, so the lipid keeps its orientation
		std::cerr << Lima::Format("Warning: could not find the head of lipid {}, using it as oriented in its .gro file\n",
			grofile.m_path.stem().string());
		return;
	}

	const Float3 lipidOrientation = (headCenter - lipidbodyCenter).norm();	

	const float alignment = lipidOrientation.dot(desiredOrientation);
	if (alignment > 0.98f) {
		return; // Already oriented correctly
	}

	Float3 rotationAxis;
	float rotationAngle;
	if (alignment < -0.98f) {
		// Upside down: the cross product vanishes, so turn half a revolution around any axis perpendicular to the lipid
		const Float3 helper = std::abs(lipidOrientation.x) < 0.9f ? Float3{ 1, 0, 0 } : Float3{ 0, 1, 0 };
		rotationAxis = lipidOrientation.cross(helper).norm();
		rotationAngle = PI;
	}
	else {
		rotationAxis = lipidOrientation.cross(desiredOrientation).norm();
		rotationAngle = std::acos(alignment);
	}
	const Float3 rotationPoint = MoleculeUtils::GeometricCenter(grofile);

	for (auto& atom : grofile.atoms) {
		atom.position -= rotationPoint;
		atom.position = Float3::rodriguesRotatation(atom.position, rotationAxis, rotationAngle);
		atom.position += rotationPoint;
	}
}




float CalculateRadius(const std::span<const Float3>& positions) {
	const Float3 mean = Statistics::Mean(positions);

	float maxDistSquared = 0;
	for (const Float3& pos : positions) {
		const float dist = (pos - mean).lenSquared();
		if (dist > maxDistSquared)
			maxDistSquared = dist;
	}
	return std::sqrt(maxDistSquared);
}

float CalculateTotalRadius(const std::vector<Float3>& positions, const std::vector<int>& partition) {
	std::vector<std::span<const Float3>> sections(partition.size());

	int head = 0;

	for (int i = 0; i < partition.size(); i++) {
		sections[i] = std::span<const Float3>(positions.data() + head, partition[i]);
		head += partition[i];
	}


	float radiusSum = 0;
	for (const auto& section : sections) {
		radiusSum += CalculateRadius(section);
	}
	return radiusSum;
}

void BruteForcePartition(const std::vector<int>& acceptableSectionSizes, const std::vector<Float3>& positions, int index, std::vector<int>& currentPartition,
	float& minRadius, std::vector<int>& bestCombination)
{
	if (index == positions.size()) {
		// Calculate radii for this partition
		double currentRadius = CalculateTotalRadius(positions, currentPartition);

		// Keep track of the minimum radius configuration
		if (currentRadius < minRadius) {
			minRadius = currentRadius;
			bestCombination = currentPartition;
		}
		return;
	}

	// Try each section size from acceptableSectionSizes
	for (int size : acceptableSectionSizes) {
		if (index + size <= positions.size()) {
			currentPartition.push_back(size);
			BruteForcePartition(acceptableSectionSizes, positions, index + size, currentPartition, minRadius, bestCombination);
			currentPartition.pop_back();
		}
	}
}

std::vector<int> DetermineAcceptableSectionsizesBasedOnAtomcount(int nAtoms) {
	int maxSize = 200;

	// SectionSizes based in 32 threads/warp
	std::vector<int> idealSectionsizes{
		64, 63, 62, 61, 60, 59, 58,
		32, 31, 30, 29, 28, 27, 26,
	};

	// Figure out if the atomcount  can be represented by these ideal sectionsizes
	{
		std::vector<bool> dp(nAtoms + 1, false);
		dp[0] = true; // base case: we can always make the sum 0

		for (int num : idealSectionsizes) {
			for (int i = nAtoms; i >= num; --i) {
				if (dp[i - num]) {
					dp[i] = true;
				}
			}
		}

		if (dp[nAtoms]) {
			return idealSectionsizes;
		}
	}

	// Otherwise, return this list that can represent all atomcounts
	{
		return std::vector<int> {
			64, 63, 62, 61, 60, 59, 58, 57, 56,
				32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21
		};
	}
}



void Lipids::OrganizeLipidIntoCompoundsizedSections(GroFile& grofile, TopologyFile::Moleculetype& molecule) {

	/*if (!topfile.GetLocalMolecules().empty()) {
		throw std::runtime_error("Cannot organise lipids with include molecules");
	}*/
	// First clear any sections previously set by LIMA 
	for (auto& atom : molecule.atoms) {
		if (atom.section_name.has_value() && atom.section_name.value() == ";lipid_section") {
			atom.section_name = std::nullopt;
		}
	}

	if (grofile.box_size == Float3{0})
		grofile.box_size = Float3{ 10.f };
	MoleculeUtils::CenterMolecule(grofile, molecule);
	OrientLipidhead(grofile);

	LimaMoleculeGraph::reorderoleculeParticlesAccoringingToSubchains(grofile, molecule);
	// Figure out how to prioritise only using the better sizes, IF that combination can be made
	std::vector<int> sectionSizes = DetermineAcceptableSectionsizesBasedOnAtomcount(grofile.atoms.size());

	std::vector<Float3> positions(grofile.atoms.size());
	for (int i = 0; i < grofile.atoms.size(); i++) {
		positions[i] = grofile.atoms[i].position;
	}

	float minRadius = FLT_MAX;
	std::vector<int> bestPartition;
	std::vector<int> currentPartition;
	BruteForcePartition(sectionSizes, positions, 0, currentPartition, minRadius, bestPartition);

	int cummulativeIndex = 0;
	for (int index : bestPartition) {
		molecule.atoms[cummulativeIndex].section_name = ";lipid_section";
		cummulativeIndex += index;
	}
}


void Lipids::_MakeLipid(const std::string& name) {
	fs::path dir = "C:/Users/Daniel/git_repo/LIMA/resources/Slipids/";
	auto groPath = dir / (name + ".gro");
	auto itpPath = dir / (name + ".itp");

	GroFile grofile{ groPath };
	TopologyFile topfile(itpPath);
	grofile.box_size = Float3{ 5.f };
	OrganizeLipidIntoCompoundsizedSections(grofile, topfile.GetMoleculeType());

	// Chicken before egg issue, cant load the lipid as it is not correctly part of the ff yet?
	{
		std::unique_ptr<Display> display = std::make_unique<Display>();
		display->Submit(0, std::make_unique<Rendering::AtomRenderTask>(grofile, true), true);
	}

	if (true) {
		grofile.printToFile();
		topfile.printToFile();
	}
}
