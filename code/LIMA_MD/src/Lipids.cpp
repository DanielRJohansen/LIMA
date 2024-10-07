#include "Lipids.h"
#include "Simulation.cuh"
#include "Statistics.h"
#include "MoleculeUtils.h"
#include "MoleculeGraph.h"
#include <format>

Lipids::Select::Select(const std::string& lipidname, const fs::path& workDir, double percentage) :
	userSupplied(fs::exists(workDir / (lipidname + ".gro")) && fs::exists(workDir / (lipidname + ".itp"))),
	lipidname(lipidname),
	percentage(percentage)
{
	const fs::path defaultLipidsDir = Filehandler::GetLimaDir() / ("resources/Slipids");

	if (!userSupplied && !fs::exists(defaultLipidsDir / (lipidname + ".itp"))) {
		throw std::runtime_error(std::format("Failed to find lipid: {}, looked here: \n\t{}\nAnd here:\n\t{}",
			lipidname, workDir.string(), defaultLipidsDir.string()));
	}

	const fs::path& lipidDir = userSupplied ? workDir : defaultLipidsDir;
	grofile = std::make_unique<GroFile>(lipidDir / (lipidname + ".gro"));
	topfile = std::make_unique<TopologyFile>(lipidDir / (lipidname + ".itp"));

	if (userSupplied) 
		OrganizeLipidIntoCompoundsizedSections(*grofile, *topfile);
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
		return; // Our estimate of the head relative to the body is too weak to orient. We wont panic but simply hope it is oriented correctly
		// TODO: We need to have a warning system in place for this for the user
	}

	const Float3 lipidOrientation = (headCenter - lipidbodyCenter).norm();	

	if (std::abs(lipidOrientation.dot(desiredOrientation)) > 0.98f) {
		return; // Already oriented correctly
	}

	const Float3 rotationAxis = lipidOrientation.cross(desiredOrientation).norm();
	const float rotationAngle = std::acos(lipidOrientation.dot(desiredOrientation));
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

#include "Display.h"

void Lipids::OrganizeLipidIntoCompoundsizedSections(GroFile& grofile, TopologyFile& topfile) {

	if (!topfile.GetLocalMolecules().empty()) {
		throw std::runtime_error("Cannot organise lipids with include molecules");
	}
	// First clear any sections previously set by LIMA 
	for (auto& atom : topfile.GetLocalAtoms()) {
		if (atom.section_name.has_value() && atom.section_name.value() == ";lipid_section") {
			atom.section_name = std::nullopt;
		}
	}

	if (grofile.box_size == Float3{0})
		grofile.box_size = Float3{ 10.f };
	MoleculeUtils::SetMoleculeCenter(grofile, grofile.box_size / 2.f);
	OrientLipidhead(grofile);

	LimaMoleculeGraph::reorderoleculeParticlesAccoringingToSubchains(grofile, topfile);
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
		topfile.GetLocalAtoms()[cummulativeIndex].section_name = ";lipid_section";
		cummulativeIndex += index;
	}
}
#include "Environment.h"

void Lipids::_MakeLipids(bool writeToFile, bool displayEachLipidAndHalt) {
	std::string path = "C:/Users/Daniel/git_repo/LIMA/resources/Slipids/";
	std::vector<std::string> targets;

	//std::unique_ptr<Display> display = displayEachLipidAndHalt ? std::make_unique<Display>(Full) : nullptr;

	for (const auto& entry : fs::directory_iterator(path)) {
		if (entry.path().extension() == ".gro") {
			std::string base_name = entry.path().stem().string();
			std::string itp_file = path + base_name + ".itp";
			if (fs::exists(itp_file)) {

				targets.push_back(base_name);
			}
		}
	}

	for (const auto& target : targets) {
		printf("Organizing %s\n", target.c_str());
		GroFile grofile{ path + target + ".gro" };
		TopologyFile topfile{ path + target + ".itp" };

		// Use the internal forcefield, so it wont matter when we end up copying the forcefield into the target dir
		assert(topfile.forcefieldIncludes.size() ==1);
		topfile.forcefieldIncludes.resize(1);
		topfile.forcefieldIncludes[0] = { TopologyFile::ForcefieldInclude{"Slipids_2020.ff/forcefield.itp"} };

		grofile.box_size = Float3{ 5.f };
	/*	if (grofile.box_size.x != grofile.box_size.y || grofile.box_size.x != grofile.box_size.z) {
			grofile.box_size = Float3{ std::max(std::max(grofile.box_size.x, grofile.box_size.y), grofile.box_size.z) };
		}*/

		OrganizeLipidIntoCompoundsizedSections(grofile, topfile);

		// Now load the lipid into a simulation. This will catch most errors we might have made in the lipid
		{
			Environment env{ grofile.m_path.parent_path(), Headless};
			SimParams params;
			params.n_steps = 2;
			params.dt = 0;
			params.data_logging_interval = 1;
			params.em_variant = true;
			env.CreateSimulation(grofile, topfile, params);
			env.run(false);
			auto sim = env.getSim();

			if (displayEachLipidAndHalt) {
				std::unique_ptr<Display> display = displayEachLipidAndHalt ? std::make_unique<Display>(Full) : nullptr; // TODO: move to top so we dont reinit every time
				display->Render(
					std::make_unique<Rendering::SimulationTask>(sim->traj_buffer->GetBufferAtStep(0), sim->box_host->compounds, sim->box_host->boxparams, 0, 0.f, ColoringMethod::GradientFromCompoundId),
					true);
			}
		}

		if (writeToFile) {
			grofile.printToFile();
			topfile.printToFile();
		}
	}
}  
