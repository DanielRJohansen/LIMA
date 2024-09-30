#include "Lipids.h"
#include "Simulation.cuh"
#include "Statistics.h"
#include "MoleculeUtils.h"
#include "MoleculeGraph.h"

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

void Lipids::_MakeLipids(std::function<void(const GroFile&, const TopologyFile&)> renderCallback, bool writeToFile) {
	std::string path = "C:/Users/Daniel/git_repo/LIMA/resources/Slipids/";
	std::vector<std::string> targets;

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
		topfile.forcefieldIncludes[0] = { TopologyFile::ForcefieldInclude{"Slipids_2020.ff/forcefield.itp", "Slipids_2020.ff/forcefield.itp"} };

		if (grofile.box_size.x != grofile.box_size.y || grofile.box_size.x != grofile.box_size.z) {
			grofile.box_size = Float3{ std::max(std::max(grofile.box_size.x, grofile.box_size.y), grofile.box_size.z) };
		}


		OrganizeLipidIntoCompoundsizedSections(grofile, topfile);

		renderCallback(grofile, topfile);
		if (writeToFile) {
			grofile.printToFile();
			topfile.printToFile();
		}
	}
}  