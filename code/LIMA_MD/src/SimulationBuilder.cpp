#include "SimulationBuilder.h"
#include "Programs.h"
#include "Format.h"

#include "BoundaryConditionPublic.h"
#include "EngineCore.h"
#include "Statistics.h"
#include "MoleculeGraph.h"
#include "MoleculeUtils.h"

#include <format>
#include <functional>
#include <fstream>
#include <algorithm>
#include <random>
#include <numeric>
#include <cfloat>
#include <numbers>




void centerMoleculeAroundOrigo(GroFile& grofile) {
	Float3 position_sum{};
	for (const auto& atom : grofile.atoms) {
		position_sum += atom.position;
	}

	const Float3 offset = position_sum / static_cast<float>(grofile.atoms.size());
	for (auto& atom : grofile.atoms) {
		atom.position -= offset;
	}
}

// Lipid structures are shared between copies of a selection, which may be building membranes concurrently,
// so the selection gets centered copies instead of centering the shared structures in place
static void UseCenteredLipidCopies(Lipids::Selection& lipidselection) {
	for (auto& lipid : lipidselection) {
		lipid.grofile = std::make_shared<GroFile>(*lipid.grofile);
		centerMoleculeAroundOrigo(*lipid.grofile);
	}
}

float constexpr MinParticlePosInDimension(const Lipids::Selection& lipidselection, int dim) {
	float minPos = FLT_MAX;
	for (const auto& lipid : lipidselection) {
		for (const auto& atom : lipid.grofile->atoms) {			
			minPos = std::min(minPos, atom.position[dim]);
		}
	}
	return minPos;
}

class RandomUniformGenerator {
	std::mt19937 generator;
	std::uniform_real_distribution<float> distribution;
public:
	RandomUniformGenerator(float min, float max, int seed = 1238971) : generator(seed), distribution(min, max) {}
	float operator()() {
		return distribution(generator);
	}
};

class RandomUniformGeneratorUnitvector {
	std::mt19937 generator;
	std::uniform_real_distribution<float> distribution;
public:
	RandomUniformGeneratorUnitvector(int seed = 1238971) : generator(seed), distribution(-1.f, 1.f) {}

	inline Float3 Generate() {
		return Float3(distribution(generator), distribution(generator), distribution(generator));
	}
	Float3 operator()() {
		Float3 val = Generate();
		while (val.lenSquared() < 1e-5)
			val = Generate();
		
		return val.norm();
	}
};


void addAtomToFile(GroFile& outputgrofile, const GroRecord& input_atom_gro, int atom_offset, int residue_offset, 
	std::function<void(Float3&)> position_transform) 
{
	outputgrofile.atoms.emplace_back(input_atom_gro);
	outputgrofile.atoms.back().gro_id += atom_offset;
	outputgrofile.atoms.back().gro_id %= 100000;	// Gro_ids only go to 99999
	outputgrofile.atoms.back().residue_number += residue_offset;
	outputgrofile.atoms.back().residue_number %= 100000; // Residue ids only go to 99999
	position_transform(outputgrofile.atoms.back().position);
}

void AddGroAndTopToGroAndTopfile(GroFile& outputgrofile, const GroFile& inputgrofile, std::function<void(Float3&)> position_transform, 
	TopologyFile& outputTopologyFile, const std::shared_ptr<TopologyFile>& inputTopology)
{
	int atomsOffset = outputgrofile.atoms.size();
	int residuenrOffset = outputgrofile.atoms.empty() ? 0 : outputgrofile.atoms.back().residue_number;

	for (const auto& atom : inputgrofile.atoms) {
		addAtomToFile(outputgrofile, atom, atomsOffset, residuenrOffset, position_transform);
	}

	if (inputTopology->GetMoleculeTypePtr() == nullptr)
		throw std::runtime_error("nullptr here");

	outputTopologyFile.AppendMoleculetype(inputTopology->GetMoleculeTypePtr(), inputTopology->forcefieldInclude);
}


void validateLipidselection(const Lipids::Selection& lipidselection) {
	double total_percentage = 0;
	for (const auto& lipid : lipidselection) {
		total_percentage += lipid.percentage;
	}
	if (std::abs(total_percentage - 100) > 0.00001f) {
		throw std::runtime_error(Lima::Format("Invalid lipid selection, did not add up to 100% {:.2f}", total_percentage));
	}

	for (const auto& lipid : lipidselection) {
		//if (lipid.grofile->atoms.size() != lipid.topfile->GetLocalAtoms().size()) {
		if (lipid.grofile->atoms.size() != lipid.topfile->GetMoleculeType().atoms.size()) {
			throw std::runtime_error(Lima::Format("BuildMembrane failed: Structure and topology file did not have the same amount of atoms. Please validate your files.\nGRO:{}\nTOP:{}",
				lipid.grofile->m_path.string(), lipid.topfile->path.string())
			);
		}
	}
}

template <typename Selection>
struct SampleSelectionRandomly {
	// Seedoffset is so we can get repeatable but different outcomes from multiple instantiations
	SampleSelectionRandomly(const Selection& selection, int seedOffset = 0)
	: selection(selection)
	, rng(34896495u + static_cast<unsigned>(seedOffset))
	{
		// Build a vector of indices reflecting the "percentage" weighting
		// (same as before)
		for (int i = 0; i < static_cast<int>(selection.size()); i++) {
			for (int j = 0; j < selection[i].percentage; j++) {
				selection_indexes.push_back(i);
			}
		}
		// Prepare a uniform distribution for picking random indices
		distribution = std::uniform_int_distribution<int>(
			0,
			static_cast<int>(selection_indexes.size()) - 1
		);
	}

	// Return a const-ref to a randomly selected lipid/atom from 'selection'
	const auto& operator()() {
		// Generate a random index to select from 'selection_indexes'
		int randomPos = distribution(rng);
		int chosenIndex = selection_indexes[randomPos];
		return selection[chosenIndex];
	}

private:
	const Selection& selection;
	std::vector<int> selection_indexes;

	std::mt19937 rng;
	std::uniform_int_distribution<int> distribution;
};

using GetNextRandomLipid = SampleSelectionRandomly<Lipids::Selection>;
using GetNextRandomParticle = SampleSelectionRandomly<AtomsSelection>;

namespace {
	constexpr float lipidDensity = 1.f / 0.59f; // [lipids/nm^2]
	constexpr int minimumInnerLeafletLipids = 12;
}

void SimulationBuilder::DistributeParticlesInBox(GroFile& grofile, TopologyFile& topfile, const AtomsSelection& particles, float minDistBetweenAnyParticle, float particlesPerNm3) 
{
	const float desiredBlockLen = 2.f;
	const int blocksPerDim = grofile.box_size.x >= 4.f
		? static_cast<int>(std::ceil(grofile.box_size.x / desiredBlockLen))
		: 1;
	const float blockLen = grofile.box_size.x / static_cast<float>(blocksPerDim);

	const float usableBlocklen = blockLen - minDistBetweenAnyParticle;

	const int particlesPerBlock = static_cast<int>(std::ceil(particlesPerNm3 * blockLen * blockLen * blockLen));


	GetNextRandomParticle getNextRandomParticle{ particles };
	RandomUniformGenerator randPos(minDistBetweenAnyParticle * 0.5f, usableBlocklen + minDistBetweenAnyParticle * 0.5f, 1238971);

	std::vector<Float3> positionsInThisBlock(particlesPerBlock);

	// Divide the box into small subblocks
	for (int z = 0; z < blocksPerDim; z++) {
		for (int y = 0; y < blocksPerDim; y++) {
			for (int x = 0; x < blocksPerDim; x++) {

				const Float3 blockStart = { x * blockLen, y * blockLen, z * blockLen };

				// For each block, create the required number of particles. Only check there's no collision inside the box, 
				// as inter-block collisions are automatically handled by the margin
				for (int relativeParticleIndex = 0; relativeParticleIndex < particlesPerBlock; ) {
					AtomtypeSelect atomtypeselect = getNextRandomParticle();

					

					const Float3 position = Float3{ randPos(), randPos(), randPos() } + blockStart;

					bool collision = false;
					for (int i = 0; i < relativeParticleIndex; i++) {
						if ((positionsInThisBlock[i] - position).len() < minDistBetweenAnyParticle) {	// No need for PBC, since it's inside blocks
							collision = true;
							break;
						}
					}

					// If no collision, add the particle to the gro and top file
					if (!collision) {
						positionsInThisBlock[relativeParticleIndex] = position;
						const int groId = grofile.atoms.empty() ? 1 : grofile.atoms.back().gro_id + 1;
						const int resNr = grofile.atoms.empty() ? 1 : grofile.atoms.back().residue_number + 1;

						grofile.atoms.emplace_back(GroRecord{ resNr, SmallString("XXX"), SmallString(atomtypeselect.atomtype.atomname), groId, position, std::nullopt });
						topfile.AppendMolecule(atomtypeselect.atomtype.atomname);

						relativeParticleIndex++;
					}
				}
			}
		}
	}
}




namespace {
struct ParticlePlaceholder {
	Float3 relPos{};	// [nm]
	bool presentInInputfile = false;	// We need to know the difference between the particles already present, and the ones we just added
	bool markedForDeletion = false;
	int waterIndex = -1;	// Which water of the equilibrated water box a new water is a copy of
};

// Cells of atoms, for finding clashes between new waters and the atoms around them. As in the engine, the box is
// split into floor(length) cells per axis, so a cell is 1 nm in an integer box and up to 2 nm otherwise; never
// shorter than the clash distances, so only the neighboring cells need searching. Periodic: a cell index at most
// one box length outside the grid wraps around
template<typename T>
class SolvationGrid {
	std::vector<std::vector<T>> nodes;
	Int3 nodesPerDim{};

public:
	const Float3 cellSize;	// [nm]

	SolvationGrid(Float3 boxSize) : nodesPerDim(BoxParams::GridSizeForLength(boxSize)),
		cellSize{ boxSize.x / nodesPerDim.x, boxSize.y / nodesPerDim.y, boxSize.z / nodesPerDim.z } {
		nodes.resize(nodesPerDim.x * nodesPerDim.y * nodesPerDim.z);
	}

	Int3 Dims() const { return nodesPerDim; }

	Float3 Origin(NodeIndex index3d) const {
		return Float3{ index3d.x * cellSize.x, index3d.y * cellSize.y, index3d.z * cellSize.z };
	}

	// The cell of a position inside the box, clamped so a position rounding to the box length stays in the last cell
	NodeIndex CellOf(Float3 position) const {
		const auto Cell = [](float p, float size, int n) { return std::clamp(static_cast<int>(std::floor(p / size)), 0, n - 1); };
		return NodeIndex{ Cell(position.x, cellSize.x, nodesPerDim.x), Cell(position.y, cellSize.y, nodesPerDim.y), Cell(position.z, cellSize.z, nodesPerDim.z) };
	}

	NodeIndex Wrap(NodeIndex index3d) const {
		index3d.x += nodesPerDim.x * ((index3d.x < 0) - (index3d.x >= nodesPerDim.x));
		index3d.y += nodesPerDim.y * ((index3d.y < 0) - (index3d.y >= nodesPerDim.y));
		index3d.z += nodesPerDim.z * ((index3d.z < 0) - (index3d.z >= nodesPerDim.z));
		return index3d;
	}

	int GetIndex(const NodeIndex& nodeindex) const {
		return nodeindex.z * nodesPerDim.x * nodesPerDim.y + nodeindex.y * nodesPerDim.x + nodeindex.x;
	}

	std::vector<T>& operator[](NodeIndex index3d) {
		return nodes[GetIndex(Wrap(index3d))];
	}
	const std::vector<T>& operator[](NodeIndex index3d) const {
		return nodes[GetIndex(Wrap(index3d))];
	}
};



void DistributeGrofileparticlesInGrid(SolvationGrid<ParticlePlaceholder>& boxgrid, const GroFile& grofile) {
	for (const auto& atom : grofile.atoms) {
		Float3 absPosHyper = atom.position;
		BoundaryConditionPublic::applyBCNM(absPosHyper, grofile.box_size, BoundaryConditionSelect::PBC);

		const NodeIndex nodeindex = boxgrid.CellOf(absPosHyper);

		// Make the positions relative to the nodeIndex
		const Float3 relPos = absPosHyper - boxgrid.Origin(nodeindex);


		boxgrid[nodeindex].emplace_back(ParticlePlaceholder{ relPos, true });
	}
}
} // namespace


// The moleculetype of a water model given as an .itp path, or as a model name looked up in LIMA's forcefields
std::shared_ptr<const TopologyFile::Moleculetype> LoadWaterMoleculetype(const std::string& waterModel, const TopologyFile& topfile) {
	fs::path itpPath;
	if (fs::path{ waterModel }.extension() == ".itp" && fs::is_regular_file(waterModel)) {
		itpPath = waterModel;
	}
	else {
		Programs::WaterModel model;
		try { model = Programs::ParseWaterModel(waterModel); }
		catch (const std::runtime_error&) {
			throw std::invalid_argument(Lima::Format("Unknown water model \"{}\": expected an existing .itp file or one of tip3p, tip4p, tips3p, tip5p, spc, spce", waterModel));
		}
		// The model's file next to the system's forcefield, which defines its atomtypes. Otherwise CHARMM27's, which has every model
		const fs::path forcefieldsDir = FileUtils::GetLimaDir() / "resources/forcefields";
		const fs::path filename{ Programs::WaterModelItpFilename(model) };
		itpPath = forcefieldsDir / "charmm27.ff" / filename;
		if (topfile.forcefieldInclude && !topfile.forcefieldInclude->filename.parent_path().empty()) {
			const fs::path candidate = forcefieldsDir / topfile.forcefieldInclude->filename.parent_path() / filename;
			if (fs::is_regular_file(candidate))
				itpPath = candidate;
		}
	}

	std::optional<TopologyFile> waterTop;
	try { waterTop.emplace(itpPath, topfile.GetDefines()); }
	catch (const std::exception& e) {
		throw std::invalid_argument(Lima::Format("Failed to load water model {}: {}", itpPath.string(), e.what()));
	}
	if (waterTop->moleculetypes.size() != 1)
		throw std::invalid_argument(Lima::Format("Water model {} must define exactly 1 moleculetype, found {}", itpPath.string(), waterTop->moleculetypes.size()));
	return waterTop->GetMoleculeTypePtr();
}

// The shape of a 3-atom water: a central atom with the 2 others bonded to it at an angle
struct WaterGeometry {
	int center = 0;
	std::array<int, 2> others{ 1, 2 };
	std::array<float, 2> bondLengths{};	// [nm]
	float angle = 0.f;					// [rad]
};

WaterGeometry GetWaterGeometry(const TopologyFile::Moleculetype& moltype) {
	if (moltype.atoms.size() != 3)
		throw std::invalid_argument(Lima::Format("Water moleculetype {} has {} atoms, solvate only supports 3-atom water (LIMA has no virtual sites)",
			moltype.name, moltype.atoms.size()));

	WaterGeometry geometry;
	if (!moltype.settles.empty()) {
		const auto& settle = moltype.settles.front();
		geometry.center = settle.oxygen;
		geometry.others = { settle.oxygen + 1, settle.oxygen + 2 };
		geometry.bondLengths = { settle.distanceOH, settle.distanceOH };
		geometry.angle = 2.f * std::asin(settle.distanceHH / (2.f * settle.distanceOH));
		return geometry;
	}

	// Flexible water: the angle's middle atom is the center, the bonds to it give the lengths
	const std::string missingGeometry = Lima::Format("Water moleculetype {} needs [ settles ], or 2 bonds and 1 angle with explicit parameters, to give its geometry", moltype.name);
	if (moltype.anglebonds.size() != 1 || !moltype.anglebonds.front().parameters || moltype.singlebonds.size() != 2)
		throw std::invalid_argument(missingGeometry);
	const auto& angle = moltype.anglebonds.front();
	geometry.center = angle.ids[1];
	geometry.others = { angle.ids[0], angle.ids[2] };
	geometry.angle = angle.parameters->theta0;
	for (int i = 0; i < 2; i++) {
		const auto bond = std::ranges::find_if(moltype.singlebonds, [&](const auto& b) {
			return (b.ids[0] == geometry.center && b.ids[1] == geometry.others[i]) || (b.ids[1] == geometry.center && b.ids[0] == geometry.others[i]); });
		if (bond == moltype.singlebonds.end() || !bond->parameters)
			throw std::invalid_argument(missingGeometry);
		geometry.bondLengths[i] = bond->parameters->b0;
	}
	return geometry;
}

void CheckWaterAtomtypesInForcefield(const TopologyFile::Moleculetype& moltype, const TopologyFile& topfile) {
	if (!topfile.forcefieldInclude)
		return;
	const auto& atomtypeEntries = topfile.forcefieldInclude->contents.GetSection(TopologySection::atomtypes);
	if (atomtypeEntries.empty())
		return;	// The forcefield's contents are unknown

	std::unordered_set<std::string> known;
	for (const auto& entry : atomtypeEntries) {
		std::istringstream stream{ entry };
		std::string name;
		stream >> name;
		known.insert(name);
	}
	std::vector<std::string> missing;
	for (const auto& atom : moltype.atoms)
		if (!known.contains(atom.type) && !std::ranges::contains(missing, atom.type))
			missing.push_back(atom.type);
	if (!missing.empty()) {
		std::string list;
		for (const auto& type : missing)
			list += (list.empty() ? "" : ", ") + type;
		throw std::invalid_argument(Lima::Format("Water moleculetype {} uses atomtypes not in the forcefield {}: {}",
			moltype.name, topfile.forcefieldInclude->filename.string(), list));
	}
}

// A cubic box of equilibrated water, tiled to solvate a system. Each water is its oxygen, wrapped into the box, and the
// directions of its 2 O-H bonds
struct EquilibratedWaterBox {
	float boxLength = 0.f;	// [nm]
	std::vector<Float3> oxygens;
	std::vector<std::array<Float3, 2>> bondDirections;
};

// GROMACS' spc216.gro: 216 SPC waters equilibrated at 300 K. Read by columns here, as GroFile would write a binary cache
// next to it, in what may be a read-only install dir
const EquilibratedWaterBox& LoadEquilibratedWaterBox() {
	static const EquilibratedWaterBox waterBox = [] {
		const fs::path path = FileUtils::GetLimaDir() / "resources/forcefields/spc216.gro";
		std::ifstream file(path);
		if (!file.is_open())
			throw std::runtime_error(Lima::Format("Failed to open {}", path.string()));

		std::string line;
		std::getline(file, line);	// Title
		std::getline(file, line);
		const int nAtoms = std::stoi(line);
		if (nAtoms % 3 != 0)
			throw std::runtime_error(Lima::Format("{} must contain 3-atom waters", path.string()));
		std::vector<Float3> positions(nAtoms);
		for (auto& position : positions) {
			std::getline(file, line);
			position = Float3{ std::stof(line.substr(20, 8)), std::stof(line.substr(28, 8)), std::stof(line.substr(36, 8)) };
		}
		std::getline(file, line);
		EquilibratedWaterBox box;
		box.boxLength = std::stof(line);

		for (int i = 0; i < nAtoms; i += 3) {
			Float3 oxygen = positions[i];
			BoundaryConditionPublic::applyBCNM(oxygen, Float3{ box.boxLength }, BoundaryConditionSelect::PBC);
			box.oxygens.push_back(oxygen);
			box.bondDirections.push_back({ (positions[i + 1] - positions[i]).norm(), (positions[i + 2] - positions[i]).norm() });
		}
		return box;
	}();
	return waterBox;
}

void SimulationBuilder::SolvateGrofile(GroFile& grofile, TopologyFile& topfile, int desiredSolventsPerNm3, const std::string& waterModel) {
	// A cell must be at least as long as the largest clash distance below, for the neighbor search to see every clash
	if (grofile.box_size.x < 1.f || grofile.box_size.y < 1.f || grofile.box_size.z < 1.f)
		throw std::invalid_argument(Lima::Format("Solvate needs a box of at least 1 nm per side, got {:.3f} x {:.3f} x {:.3f} nm",
			grofile.box_size.x, grofile.box_size.y, grofile.box_size.z));

	// Solvating an empty box: the topology gets the forcefield of the water model, and a system to add the water to
	if (!topfile.forcefieldInclude) {
		topfile.forcefieldInclude = TopologyFile{ FileUtils::GetLimaDir() / "resources/forcefields/Slipids_2020.ff/forcefield.itp" }.forcefieldInclude;
		topfile.forcefieldInclude->filename = "Slipids_2020.ff/forcefield.itp";
	}
	const auto waterMoltype = LoadWaterMoleculetype(waterModel, topfile);
	const WaterGeometry waterGeometry = GetWaterGeometry(*waterMoltype);
	CheckWaterAtomtypesInForcefield(*waterMoltype, topfile);
	if (topfile.moleculetypes.contains(waterMoltype->name) && topfile.moleculetypes.at(waterMoltype->name)->atoms != waterMoltype->atoms)
		throw std::invalid_argument(Lima::Format("The topology already has a different moleculetype named {}", waterMoltype->name));

	const int nAtomsInput = grofile.atoms.size();
	const int startResidueId = grofile.atoms.empty() ? 1 : grofile.atoms.back().residue_number + 1;

	SolvationGrid<ParticlePlaceholder> boxgrid{ grofile.box_size };
	const Int3 gridDim = boxgrid.Dims();


	//SolvationGrid<Float3> nonSolventPositions{ gridDim };	// This is used to store the positions of the particles that are not solvents, so we can remove them later
	//for (auto elem : grofile.atoms) {				
	//	NodeIndex nodeindex = NodeIndex{ static_cast<int>(std::floor(elem.position.x)), static_cast<int>(std::floor(elem.position.y)), static_cast<int>(std::floor(elem.position.z)) };
	//	Float3 relpos = elem.position - Float3{ static_cast<float>(nodeindex.x), static_cast<float>(nodeindex.y), static_cast<float>(nodeindex.z) };

	//	nonSolventPositions[nodeindex].emplace_back(relpos);
	//}

	DistributeGrofileparticlesInGrid(boxgrid, grofile);

	// Minimum distance from a new water to another new water. Below the closest O-O distance in the equilibrated water
	// box (0.251 nm), so it only removes waters that clash where the tiles meet and across the box's periodic boundary
	const float distanceThreshold = 0.23f;	// [nm]
	// Minimum distance from a new water to an atom of the input system. Larger than between waters, so waters
	// are not placed in the gaps between packed molecules, eg. in the hydrophobic core of a membrane
	const float soluteDistanceThreshold = 0.35f;	// [nm]

	// Tile the box with copies of the equilibrated water box, scaled to the desired density, keeping the waters whose
	// oxygen is inside the box. The waters keep their orientations, so they start out hydrogen bonded
	const EquilibratedWaterBox& waterBox = LoadEquilibratedWaterBox();
	const float nativeDensity = waterBox.oxygens.size() / (waterBox.boxLength * waterBox.boxLength * waterBox.boxLength);	// [1/nm^3]
	const float scale = std::cbrt(nativeDensity / static_cast<float>(desiredSolventsPerNm3));
	if (0.251f * scale < distanceThreshold)
		throw std::invalid_argument(Lima::Format("Solvent density {}/nm^3 is too high, the maximum is {}", desiredSolventsPerNm3,
			static_cast<int>(nativeDensity * std::pow(0.251f / distanceThreshold, 3.f))));
	const float tileLength = waterBox.boxLength * scale;
	const Int3 tileDim{ static_cast<int>(std::ceil(grofile.box_size.x / tileLength)), static_cast<int>(std::ceil(grofile.box_size.y / tileLength)),
		static_cast<int>(std::ceil(grofile.box_size.z / tileLength)) };

	for (int x = 0; x < tileDim.x; x++) {
		for (int y = 0; y < tileDim.y; y++) {
			for (int z = 0; z < tileDim.z; z++) {
				const Float3 tileOffset = Float3{ static_cast<float>(x), static_cast<float>(y), static_cast<float>(z) } * tileLength;
				for (int i = 0; i < static_cast<int>(waterBox.oxygens.size()); i++) {
					const Float3 position = waterBox.oxygens[i] * scale + tileOffset;
					if (position.x >= grofile.box_size.x || position.y >= grofile.box_size.y || position.z >= grofile.box_size.z)
						continue;
					const NodeIndex nodeindex = boxgrid.CellOf(position);
					const Float3 relPos = position - boxgrid.Origin(nodeindex);
					boxgrid[nodeindex].emplace_back(ParticlePlaceholder{ relPos, false, false, i });
				}
			}
		}
	}

	// Now mark all particles too close to another for deletion, if said particle is the "lower" id/block compared to the other
	for (int x = 0; x < gridDim.x; x++) {
		for (int y = 0; y < gridDim.y; y++) {
			for (int z = 0; z < gridDim.z; z++) {
				const NodeIndex nodeindex = NodeIndex{ x, y, z };
				auto& particles = boxgrid[nodeindex];
				
				// First search through all solvent particles in this block
				for (int pid = 0; pid < particles.size(); pid++) {
					if (particles[pid].presentInInputfile)	// Cant delete any particles we did not place
						continue;

					for (int otherId = 0; otherId < particles.size(); otherId++) {
						if (pid == otherId)
							continue;
						if (particles[otherId].markedForDeletion)
							continue;

						if (otherId < pid && !particles[otherId].presentInInputfile)
							continue;


						const float threshold = particles[otherId].presentInInputfile ? soluteDistanceThreshold : distanceThreshold;
						if ((particles[pid].relPos - particles[otherId].relPos).len() < threshold) {
							particles[pid].markedForDeletion = true;
							break;
						}
					}
				}
				// Now search nonsolvents in block



				// Now search through all surrounding blocks
				for (int offsetX = -1; offsetX < 2; offsetX++) {
					for (int offsetY = -1; offsetY < 2; offsetY++) {
						for (int offsetZ = -1; offsetZ < 2; offsetZ++) {
							if (offsetX == 0 && offsetY == 0 && offsetZ == 0)
								continue;

							const NodeIndex otherNodeIndex = boxgrid.Wrap(NodeIndex{ x + offsetX, y + offsetY, z + offsetZ });

							const Float3 relPosOffsetOther = boxgrid.Origin(NodeIndex{ offsetX, offsetY, offsetZ });

							for (const auto& queryParticle : boxgrid[otherNodeIndex]) {
								if (queryParticle.markedForDeletion)
									continue;
								const Float3 queryParticlePosition = queryParticle.relPos + relPosOffsetOther;

								for (auto& particle : particles) {
									if (particle.presentInInputfile || particle.markedForDeletion)
										continue;

									// If our particle is of the greater nodeindes, then we wont remove it, so we can continue
									if (!queryParticle.presentInInputfile && boxgrid.GetIndex(otherNodeIndex) < boxgrid.GetIndex(nodeindex))
										continue;

									const float threshold = queryParticle.presentInInputfile ? soluteDistanceThreshold : distanceThreshold;
									if ((particle.relPos - queryParticlePosition).len() < threshold) {
										particle.markedForDeletion = true;
									}
								}
							}

						}
					}
				}



			}
		}
	}
	
	// The outer atoms are placed in the plane and around the bisector of the equilibrated water's O-H bonds, at the
	// bond lengths and angle of the water model
	const float cosHalfAngle = std::cos(waterGeometry.angle * 0.5f);
	const float sinHalfAngle = std::sin(waterGeometry.angle * 0.5f);

	int atomCount = nAtomsInput;
	int solventCount = 0;
	for (int x = 0; x < gridDim.x; x++) {
		for (int y = 0; y < gridDim.y; y++) {
			for (int z = 0; z < gridDim.z; z++) {
				const NodeIndex nodeindex = NodeIndex{ x, y, z };
				for (const auto& solvent : boxgrid[nodeindex]) {
					if (solvent.markedForDeletion || solvent.presentInInputfile)
						continue;


					const auto& bonds = waterBox.bondDirections[solvent.waterIndex];
					const Float3 bisector = (bonds[0] + bonds[1]).norm();
					const Float3 inPlane = (bonds[0] - bonds[1]).norm();
					const Float3 h1Pos = (bisector * cosHalfAngle + inPlane * sinHalfAngle) * waterGeometry.bondLengths[0];
					const Float3 h2Pos = (bisector * cosHalfAngle - inPlane * sinHalfAngle) * waterGeometry.bondLengths[1];

					const Float3 blockOffset = boxgrid.Origin(nodeindex);
					std::array<Float3, 3> positions;
					positions[waterGeometry.center] = solvent.relPos + blockOffset;
					positions[waterGeometry.others[0]] = solvent.relPos + blockOffset + h1Pos;
					positions[waterGeometry.others[1]] = solvent.relPos + blockOffset + h2Pos;
					for (int i = 0; i < 3; i++) {
						const auto& atom = waterMoltype->atoms[i];
						grofile.atoms.push_back(GroRecord{ (solventCount + startResidueId) % 100000, SmallString(atom.residue), SmallString(atom.atomname),
							(atomCount + i + 1) % 100000, positions[i], std::nullopt });
					}

					atomCount += 3;
					solventCount++;
				}
			}
		}
	}

	if (!topfile.HasSystem())
		topfile.SetSystem("Water");
	topfile.AppendMoleculetype(waterMoltype);
	// AppendMoleculetype also appends the first molecule.
	for (int i = 1; i < solventCount; i++) {
		topfile.AppendMolecule(waterMoltype->name);
	}
}

void SimulationBuilder::InsertSubmoleculeInSimulation(GroFile& targetGrofile, TopologyFile& targetTopol,
	GroFile& submolGro, const std::shared_ptr<TopologyFile>& submolTop, Float3 targetCenter)
{
	MoleculeUtils::CenterMolecule(submolGro, submolTop->GetMoleculeType());
	const Float3 molCenter = MoleculeUtils::GeometricCenter(submolGro);

	std::function<void(Float3&)> position_transform = [&](Float3& pos) {
		pos -= molCenter;
		pos += targetCenter;
	};

	AddGroAndTopToGroAndTopfile(targetGrofile, submolGro, position_transform, targetTopol, submolTop);
}

void SimulationBuilder::InsertSystem(GroFile& targetGrofile, TopologyFile& targetTopol, const GroFile& sourceGrofile,
	const TopologyFile& sourceTopol, Float3 targetCenter,
	const std::function<bool(const TopologyFile::Moleculetype&)>& keep)
{
	// The source atoms of the kept molecules, in topology order
	struct Range { size_t begin, end; };
	std::vector<Range> kept;
	std::vector<const TopologyFile::MoleculeEntry*> keptEntries;
	size_t atom = 0;
	Float3 centerSum{};
	size_t nKeptAtoms = 0;
	for (const auto& entry : sourceTopol.GetSystem().molecules) {
		const size_t nAtoms = entry.moleculetype->atoms.size() * entry.count;
		if (atom + nAtoms > sourceGrofile.atoms.size())
			throw std::runtime_error("InsertSystem: the source topology has more atoms than its coordinates");
		if (keep(*entry.moleculetype)) {
			kept.push_back({ atom, atom + nAtoms });
			keptEntries.push_back(&entry);
			for (size_t i = atom; i < atom + nAtoms; i++) centerSum += sourceGrofile.atoms[i].position;
			nKeptAtoms += nAtoms;
		}
		atom += nAtoms;
	}
	if (nKeptAtoms == 0)
		throw std::invalid_argument("InsertSystem: no molecules to insert");
	const Float3 offset = targetCenter - centerSum / static_cast<float>(nKeptAtoms);

	const int atomsOffset = static_cast<int>(targetGrofile.atoms.size());
	const int residueOffset = targetGrofile.atoms.empty() ? 0 : targetGrofile.atoms.back().residue_number;
	for (const Range& range : kept)
		for (size_t i = range.begin; i < range.end; i++)
			addAtomToFile(targetGrofile, sourceGrofile.atoms[i], atomsOffset, residueOffset, [&](Float3& position) { position += offset; });

	for (const auto* entry : keptEntries) {
		if (const auto existing = targetTopol.moleculetypes.find(entry->name); existing == targetTopol.moleculetypes.end())
			targetTopol.moleculetypes.insert({ entry->name, std::make_shared<TopologyFile::Moleculetype>(*entry->moleculetype) });
		else if (existing->second->atoms.size() != entry->moleculetype->atoms.size())
			throw std::runtime_error(Lima::Format("InsertSystem: both systems have a different moleculetype named {}", entry->name));
		for (int i = 0; i < entry->count; i++)
			targetTopol.AppendMolecule(entry->name);
	}

	if (sourceTopol.forcefieldInclude) {
		if (!targetTopol.forcefieldInclude)
			targetTopol.forcefieldInclude = sourceTopol.forcefieldInclude;
		else {
			for (const TopologySection section : { atomtypes, pairtypes, bondtypes, constainttypes, angletypes, dihedraltypes, impropertypes, cmaptypes })
				for (const std::string& line : sourceTopol.forcefieldInclude->contents.GetSection(section))
					targetTopol.forcefieldInclude->AddEntry(section, line);
			// No longer a stock forcefield, so printing the topology must write its contents rather than include the stock file
			targetTopol.forcefieldInclude->filename = "forcefield.itp";
		}
	}
}

void SimulationBuilder::InsertSubmoleculesInSimulation(GroFile& targetGrofile, TopologyFile& targetTopol,
	GroFile& submolGro, const std::shared_ptr<TopologyFile>& submolTop, int nMoleculesToInsert, bool rotateRandomly) 
{
	MoleculeUtils::CenterMolecule(submolGro, submolTop->GetMoleculeType());

	if (submolTop->moleculetypes.size() > 1)
		throw std::invalid_argument("Source topology contains more than 1 moleculetype, which is not allowed. Check your file for possible #include with moleculetype defidisnitions");

	const Float3 molCenter = MoleculeUtils::GeometricCenter(submolGro);
	const float molRadius = MoleculeUtils::Radius(submolGro, molCenter) * 1.1f;

	RandomUniformGenerator genRandomAngle(-PI, PI);	
	RandomUniformGenerator genTranslation(molRadius, targetGrofile.box_size.x - molRadius);
	assert(targetGrofile.box_size.x == targetGrofile.box_size.y && targetGrofile.box_size.x == targetGrofile.box_size.z);

	targetGrofile.atoms.reserve(targetGrofile.atoms.size() + nMoleculesToInsert * submolGro.atoms.size());

	for (int i = 0; i < nMoleculesToInsert; i++) {
		Float3 randomTranslation = Float3{ genTranslation(), genTranslation() , genTranslation() };
		Float3 randomRotation = Float3{genRandomAngle(), genRandomAngle(), genRandomAngle()};

		std::function<void(Float3&)> position_transform = [=](Float3& pos) {
			pos -= molCenter;
			if (rotateRandomly) {
				pos = Float3::rodriguesRotatation(pos, Float3(0,0,1), randomRotation.z);
				pos = Float3::rodriguesRotatation(pos, Float3(0,1,0), randomRotation.y);
				pos = Float3::rodriguesRotatation(pos, Float3(1,0,0), randomRotation.x);
			}
			
			pos += randomTranslation;
			};

		AddGroAndTopToGroAndTopfile(targetGrofile, submolGro, position_transform, targetTopol, submolTop);
	}
}

void SimulationBuilder::InsertSubmoleculesOnSphere(
	GroFile& targetGrofile,
	TopologyFile& targetTopol,
	Lipids::Selection lipidselection,
	int nMoleculesToInsert,
	float sphereRadius,
	const Float3& sphereCenter
)
{
	RandomUniformGenerator genRandomAngle(-PI, PI);

	UseCenteredLipidCopies(lipidselection);


	GetNextRandomLipid genNextRandomLipid{ lipidselection };
	RandomUniformGenerator genTranslation(-0.5f, 0.5f);

	// Use Fibonacci lattice to evenly distribute points on a sphere
	const float phi = (1.0f + std::sqrt(5.0f)) / 2.0f; // Golden ratio

	for (int i = 0; i < nMoleculesToInsert; i++) {
		float z = 1.0f - (2.0f * i) / static_cast<float>(std::max(nMoleculesToInsert,2) - 1); // z-coordinate
		float radius = std::sqrt(1.0f - z * z); // radius for current z slice

		float theta = 2.0f * PI * i / phi; // angle in xy-plane

		Float3 translationToPointOnSphere = Float3{
			sphereRadius * radius * cos(theta),
			sphereRadius * radius * sin(theta),
			sphereRadius * z
		} + sphereCenter;

		// Calculate the outward normal vector at this point on the sphere
		Float3 outwardNormal = Float3{
			radius * cos(theta),
			radius * sin(theta),
			z
		};
		outwardNormal.norm(); // Ensure the normal vector is a unit vector

		// Determine the rotation needed to align the molecule's up direction (0, 0, 1) with the outward normal
		Float3 currentUp = Float3{ 0.0f, 0.0f, 1.0f };
		
		Float3 rotationAxis;
		float rotationAngle;

		if (std::abs(outwardNormal.z - (-1.0f)) < 1e-6) { // Directly downward (antiparallel case)
			// Rotate 180 degrees around x-axis or y-axis if vectors are opposite
			rotationAxis = Float3{ 1.0f, 0.0f, 0.0f };
			rotationAngle = PI;
		}
		else {
			rotationAxis = currentUp.cross(outwardNormal);
			if (rotationAxis.len() < 1e-3) {
				rotationAxis = Float3{ 1.0f, 0.0f, 0.0f }; // Any perpendicular vector if they're parallel
			}

			rotationAxis = rotationAxis.norm();
			rotationAngle = std::acos(currentUp.dot(outwardNormal)); // Angle between current up and outward normal
		}
		//Float3 rotationNoise = Float3{ genRandomAngle(), genRandomAngle(), genRandomAngle() } / 2.f;
		//rotationAxis = (rotationAxis + rotationNoise).norm();

		Float3 translationNoise = Float3{ genTranslation(), genTranslation(), genTranslation() };



		const Lipids::Select& lipid = genNextRandomLipid();

		const float selfRotAngle = genRandomAngle();
		std::function<void(Float3&)> position_transform = [&](Float3& pos) {
			pos = Float3::rodriguesRotatation(pos, Float3{0,0,1}, selfRotAngle); // Randomly rotate the molecule along its own axis

			pos = Float3::rodriguesRotatation(pos, rotationAxis, rotationAngle); // Rotate around the calculated axis by the angl

			pos += translationToPointOnSphere;
			
			pos += translationNoise;
			};

		AddGroAndTopToGroAndTopfile(targetGrofile, *lipid.grofile, position_transform, targetTopol, lipid.topfile);
	}
}








struct QueuedInsertion {
	GroFile& grofile;
	std::function<void(Float3&)> positionTransform;
	std::shared_ptr<TopologyFile> topfile;
};

static std::pair<Float3, float> RotationFromZAxisTo(const Float3& targetDirection);

// Thermal membrane roughness is spatially correlated. Independent large
// displacements would create holes and unrealistic bilayer-thickness changes.
// These periodic low-frequency modes give a nanoscale patch sub-nanometre
// undulations while remaining continuous across periodic box boundaries.
static float PlanarMidplaneDisplacement(float x, float y, const Float3& boxSize) {
	const float shortestSide = std::min(boxSize.x, boxSize.y);
	const float amplitude = std::clamp(shortestSide * 0.025f, 0.10f, 0.7f); // [nm]
	const float xPhase = 2.f * std::numbers::pi_v<float> * x / boxSize.x;
	const float yPhase = 2.f * std::numbers::pi_v<float> * y / boxSize.y;
	return amplitude * (
		0.55f * std::sin(xPhase + 0.37f) * std::cos(yPhase - 0.61f)
		+ 0.30f * std::cos(2.f * xPhase - 0.83f) * std::sin(yPhase + 0.29f)
		+ 0.15f * std::sin(xPhase + yPhase + 1.17f));
}

static Float3 PlanarSurfaceNormal(float x, float y, const Float3& boxSize) {
	constexpr float sampleDistance = 0.01f; // [nm]
	const float dhdx = (PlanarMidplaneDisplacement(x + sampleDistance, y, boxSize)
		- PlanarMidplaneDisplacement(x - sampleDistance, y, boxSize)) / (2.f * sampleDistance);
	const float dhdy = (PlanarMidplaneDisplacement(x, y + sampleDistance, boxSize)
		- PlanarMidplaneDisplacement(x, y - sampleDistance, boxSize)) / (2.f * sampleDistance);
	return Float3{ -dhdx, -dhdy, 1.f }.norm();
}

static void CreatePlanarMembrane(GroFile& grofile, TopologyFile& topfile,
	const Lipids::Selection& lipidselection, float membraneCenter, int randomSeed) {

	const float lowestZpos = MinParticlePosInDimension(lipidselection, 2);
	const float n_lipids_total = lipidDensity * grofile.box_size.x * grofile.box_size.y; // (per side)
	const int lipidsPerDimx = static_cast<int>(std::ceil(sqrtf(n_lipids_total)));
	const int lipidsPerDimy = static_cast<int>(std::ceil(n_lipids_total / static_cast<float>(lipidsPerDimx)));

	const float distPerX = grofile.box_size.x / static_cast<float>(lipidsPerDimx);
	const float distPerY = grofile.box_size.y / static_cast<float>(lipidsPerDimy);


	const float interLipidLayerSpaceHalf = 0.01f; // [nm]

	RandomUniformGenerator genRandomAngle(-PI, PI, 1238971 + randomSeed);
	GetNextRandomLipid getNextRandomLipid{ lipidselection, randomSeed };
	// Small independent protrusions are layered on top of the shared mid-plane
	// undulation. Keeping these below 0.04 nm avoids tearing the initial bilayer.
	RandomUniformGenerator topLeafletProtrusion(-0.035f, 0.035f, 571923 + randomSeed);
	RandomUniformGenerator bottomLeafletProtrusion(-0.035f, 0.035f, 927531 + randomSeed);

	std::map<std::string, std::vector<QueuedInsertion>> queuedInsertions; // Must be ordered, so we get the same sequence each time
	for (auto& lipid : lipidselection) {
		queuedInsertions[lipid.lipidname] = std::vector<QueuedInsertion>();
	}

	int nLipidsInserted = 0;
	for (int x = 0; x < lipidsPerDimx; x++) {
		const float packingOffset = x % 2 == 0 ? 0.f : distPerX / 2.f;

		for (int y = 0; y < lipidsPerDimy; y++) {

			if (nLipidsInserted == static_cast<int>(n_lipids_total))
				break;


			const float centerX = static_cast<float>(x) * distPerX + distPerX / 2.f;
			const float centerY = static_cast<float>(y) * distPerY + distPerY / 4.f + packingOffset;
			const float midplaneDisplacement = PlanarMidplaneDisplacement(centerX, centerY, grofile.box_size);
			const Float3 surfaceNormal = PlanarSurfaceNormal(centerX, centerY, grofile.box_size);

			// Insert top lipid
			{
				const Lipids::Select& inputlipid = getNextRandomLipid();

				const Float3 lipidCenter = Float3{
					centerX,
					centerY,
					membraneCenter + midplaneDisplacement + std::abs(lowestZpos)
						+ interLipidLayerSpaceHalf + topLeafletProtrusion()
				};

				const float randomRot = genRandomAngle();
				const auto [normalRotationAxis, normalRotationAngle] = RotationFromZAxisTo(surfaceNormal);
				std::function<void(Float3&)> position_transform =
					[lipidCenter, randomRot, normalRotationAxis, normalRotationAngle](Float3& pos) {
					pos = Float3::rodriguesRotatation(pos, Float3{ 0,0,1 }, randomRot);
					pos = Float3::rodriguesRotatation(pos, normalRotationAxis, normalRotationAngle);
					pos += lipidCenter;
					};

				queuedInsertions.at(inputlipid.lipidname).emplace_back(QueuedInsertion{ *inputlipid.grofile, position_transform, inputlipid.topfile });
			}

			// Insert bottom lipid
			{
				const Lipids::Select& inputlipid = getNextRandomLipid();

				const Float3 lipidCenter = Float3{
					centerX,
					centerY,
					membraneCenter + midplaneDisplacement - std::abs(lowestZpos)
						- interLipidLayerSpaceHalf + bottomLeafletProtrusion()
				};

				const float randomRot = genRandomAngle();
				const auto [normalRotationAxis, normalRotationAngle] = RotationFromZAxisTo(-surfaceNormal);
				std::function<void(Float3&)> position_transform =
					[lipidCenter, randomRot, normalRotationAxis, normalRotationAngle](Float3& pos) {
					pos = Float3::rodriguesRotatation(pos, Float3{ 0,0,1 }, randomRot);
					pos = Float3::rodriguesRotatation(pos, normalRotationAxis, normalRotationAngle);
					pos += lipidCenter;
					};

				//AddGroAndTopToGroAndTopfile(grofile, *inputlipid.grofile, position_transform,
				//	topfile, inputlipid.topfile);
				queuedInsertions.at(inputlipid.lipidname).emplace_back(QueuedInsertion{ *inputlipid.grofile, position_transform, inputlipid.topfile });
			}

			nLipidsInserted++;
		}
	}

// This doesnt compile with GCC, but it enables execution:par, no?
/*	const int totalIncoming = std::reduce(queuedInsertions.begin(), queuedInsertions.end(), 0, [](int sum, const auto& pair) {
		const int atomsPerLipid = pair.second.front().grofile.atoms.size();
		const int nLipids = pair.second.size();
		return sum + nLipids*atomsPerLipid;
		});*/
	const int totalIncoming = std::accumulate(queuedInsertions.begin(), queuedInsertions.end(), 0, [](int sum, const auto& pair) {
		if (pair.second.empty())
			return sum;
		const int atomsPerLipid = pair.second.front().grofile.atoms.size();
		const int nLipids = pair.second.size();
    return sum + nLipids * atomsPerLipid;
});

	grofile.atoms.reserve(grofile.atoms.size() + totalIncoming);

	for (const auto& [_, lipidType] : queuedInsertions) {
		for (const QueuedInsertion& queuedElement : lipidType) {
			AddGroAndTopToGroAndTopfile(grofile, queuedElement.grofile, queuedElement.positionTransform, topfile, queuedElement.topfile);
		}
	}

}

float SimulationBuilder::MinimumSphereRadius(const Lipids::Selection& lipidselection) {
	if (lipidselection.empty())
		throw std::invalid_argument("Cannot determine a membrane sphere radius without lipids.");

	float leafletHalfThickness = 0.f;
	for (const auto& lipid : lipidselection) {
		if (lipid.grofile->atoms.empty())
			throw std::invalid_argument("Cannot build a membrane from an empty lipid structure.");

		Float3 center{};
		for (const auto& atom : lipid.grofile->atoms)
			center += atom.position;
		center = center / static_cast<float>(lipid.grofile->atoms.size());
		for (const auto& atom : lipid.grofile->atoms)
			leafletHalfThickness = std::max(leafletHalfThickness, center.z - atom.position.z);
	}

	const float minimumInnerRadius = std::sqrt(
		static_cast<float>(minimumInnerLeafletLipids) /
		(4.f * std::numbers::pi_v<float> * lipidDensity));
	return leafletHalfThickness + minimumInnerRadius;
}

static std::pair<Float3, float> RotationFromZAxisTo(const Float3& targetDirection) {
	const Float3 zAxis{ 0.f, 0.f, 1.f };
	const Float3 target = targetDirection.norm();
	const float cosine = std::clamp(zAxis.dot(target), -1.f, 1.f);
	if (cosine < -1.f + 1e-6f)
		return { Float3{ 1.f, 0.f, 0.f }, PI };

	Float3 axis = zAxis.cross(target);
	if (axis.lenSquared() < 1e-8f)
		axis = Float3{ 1.f, 0.f, 0.f };
	else
		axis = axis.norm();
	return { axis, std::acos(cosine) };
}

// A sum of low-order, zero-mean spherical modes approximates the correlated
// thermal shape fluctuations of a vesicle without changing its mean radius.
static float ClosedSurfaceDisplacement(const Float3& direction, float characteristicRadius) {
	const float amplitude = std::clamp(characteristicRadius * 0.035f, 0.10f, 0.60f); // [nm]
	const float quadrupole = 0.5f * (3.f * direction.z * direction.z - 1.f);
	return amplitude * (
		0.50f * quadrupole
		+ 0.28f * (direction.x * direction.x - direction.y * direction.y)
		+ 0.14f * (2.f * direction.x * direction.y)
		+ 0.08f * (2.f * direction.y * direction.z));
}

static Float3 FibonacciDirection(float index, int count) {
	const float goldenAngle = std::numbers::pi_v<float> * (3.f - std::sqrt(5.f));
	const float z = 1.f - 2.f * (index + 0.5f) / static_cast<float>(count);
	const float xyRadius = std::sqrt(std::max(0.f, 1.f - z * z));
	const float theta = goldenAngle * index;
	return Float3{ xyRadius * std::cos(theta), xyRadius * std::sin(theta), z };
}

// Surface-area Jacobian for the parameterization
// (a sqrt(1-u^2) cos(theta), b sqrt(1-u^2) sin(theta), c u).
static double EllipsoidSurfaceJacobian(double u, double theta, const Float3& radii) {
	const double a = radii.x;
	const double b = radii.y;
	const double c = radii.z;
	const double oneMinusUSquared = std::max(0.0, 1.0 - u * u);
	const double cosine = std::cos(theta);
	const double sine = std::sin(theta);
	return std::sqrt(
		b * b * c * c * oneMinusUSquared * cosine * cosine
		+ a * a * c * c * oneMinusUSquared * sine * sine
		+ a * a * b * b * u * u);
}

// Generate a deterministic equal-area point set. The old implementation
// resampled a continuous Fibonacci spiral at fractional indices. Even a small
// departure from a sphere then produced adjacent angles much closer than the
// golden angle, creating bands and sometimes overlapping lipids.
//
// Here each integer Fibonacci index selects one equal-area meridional band and
// one golden-angle azimuthal quantile. Both quantiles are transformed through
// the ellipsoid's surface-area Jacobian, so triaxial ellipsoids remain uniform
// without sacrificing the separation properties of the integer lattice.
static std::vector<Float3> EllipsoidSurfaceDirections(int count, const Float3& radii) {
	const float longestAxis = std::max({ radii.x, radii.y, radii.z });
	const float shortestAxis = std::min({ radii.x, radii.y, radii.z });
	if (longestAxis - shortestAxis <= longestAxis * 1e-6f) {
		std::vector<Float3> directions;
		directions.reserve(count);
		for (int i = 0; i < count; ++i)
			directions.push_back(FibonacciDirection(static_cast<float>(i), count));
		return directions;
	}

	constexpr int meridionalIntervals = 2048;
	constexpr int azimuthalIntervals = 256;
	constexpr double twoPi = 2.0 * std::numbers::pi_v<double>;
	const double goldenAngle = std::numbers::pi_v<double> * (3.0 - std::sqrt(5.0));

	// Integrate out theta to obtain the marginal surface-area CDF for u.
	std::vector<double> meridionalCdf(meridionalIntervals + 1, 0.0);
	auto marginalAreaDensity = [&radii](double u) {
		double sum = 0.0;
		for (int j = 0; j < azimuthalIntervals; ++j) {
			const double theta = twoPi * (static_cast<double>(j) + 0.5)
				/ static_cast<double>(azimuthalIntervals);
			sum += EllipsoidSurfaceJacobian(u, theta, radii);
		}
		return sum / static_cast<double>(azimuthalIntervals);
	};

	double previousDensity = marginalAreaDensity(-1.0);
	for (int j = 1; j <= meridionalIntervals; ++j) {
		const double u = -1.0 + 2.0 * static_cast<double>(j)
			/ static_cast<double>(meridionalIntervals);
		const double density = marginalAreaDensity(u);
		meridionalCdf[j] = meridionalCdf[j - 1] + previousDensity + density;
		previousDensity = density;
	}
	const double totalMeridionalWeight = meridionalCdf.back();

	std::vector<Float3> directions;
	directions.reserve(count);
	for (int i = 0; i < count; ++i) {
		// Descending u matches the pole ordering of FibonacciDirection.
		const double meridionalTarget = totalMeridionalWeight
			* (1.0 - (static_cast<double>(i) + 0.5) / static_cast<double>(count));
		const auto upper = std::lower_bound(
			meridionalCdf.begin() + 1, meridionalCdf.end(), meridionalTarget);
		const int upperIndex = static_cast<int>(upper - meridionalCdf.begin());
		const int lowerIndex = upperIndex - 1;
		const double intervalWeight = meridionalCdf[upperIndex] - meridionalCdf[lowerIndex];
		const double fraction = intervalWeight > 0.0
			? (meridionalTarget - meridionalCdf[lowerIndex]) / intervalWeight : 0.5;
		const double u = -1.0 + 2.0
			* (static_cast<double>(lowerIndex) + fraction)
			/ static_cast<double>(meridionalIntervals);

		// Convert the integer golden-angle phase to the conditional equal-area
		// azimuth at this u. Midpoint quadrature avoids bias at bin boundaries.
		double azimuthalWeight = 0.0;
		for (int j = 0; j < azimuthalIntervals; ++j) {
			const double theta = twoPi * (static_cast<double>(j) + 0.5)
				/ static_cast<double>(azimuthalIntervals);
			azimuthalWeight += EllipsoidSurfaceJacobian(u, theta, radii);
		}
		const double phase = std::fmod(goldenAngle * static_cast<double>(i), twoPi) / twoPi;
		const double azimuthalTarget = phase * azimuthalWeight;
		double cumulativeWeight = 0.0;
		double theta = 0.0;
		for (int j = 0; j < azimuthalIntervals; ++j) {
			const double midpoint = twoPi * (static_cast<double>(j) + 0.5)
				/ static_cast<double>(azimuthalIntervals);
			const double weight = EllipsoidSurfaceJacobian(u, midpoint, radii);
			if (cumulativeWeight + weight >= azimuthalTarget) {
				const double withinBin = weight > 0.0
					? (azimuthalTarget - cumulativeWeight) / weight : 0.5;
				theta = twoPi * (static_cast<double>(j) + withinBin)
					/ static_cast<double>(azimuthalIntervals);
				break;
			}
			cumulativeWeight += weight;
		}

		const double xyRadius = std::sqrt(std::max(0.0, 1.0 - u * u));
		directions.emplace_back(
			static_cast<float>(xyRadius * std::cos(theta)),
			static_cast<float>(xyRadius * std::sin(theta)),
			static_cast<float>(u));
	}
	return directions;
}

static Float3 EllipsoidNormal(const Float3& direction, const Float3& radii) {
	return Float3{
		direction.x / radii.x,
		direction.y / radii.y,
		direction.z / radii.z
	}.norm();
}

static float EllipsoidSurfaceArea(const Float3& radii) {
	// Knud Thomsen's approximation is well below the uncertainty in the target
	// lipid area and is exact for spheres.
	constexpr float p = 1.6075f;
	const float mean = (
		std::pow(radii.x * radii.y, p)
		+ std::pow(radii.x * radii.z, p)
		+ std::pow(radii.y * radii.z, p)) / 3.f;
	return 4.f * std::numbers::pi_v<float> * std::pow(mean, 1.f / p);
}

static void QueueEllipsoidalLeaflet(
	std::map<std::string, std::vector<QueuedInsertion>>& queuedInsertions,
	const Lipids::Selection& lipidselection,
	const Float3& center,
	const Float3& midSurfaceRadii,
	float normalOffset,
	int lipidCount,
	bool pointOutwards,
	int randomSeedOffset)
{
	GetNextRandomLipid getNextRandomLipid{ lipidselection, randomSeedOffset };
	RandomUniformGenerator randomRotation(-PI, PI, 1238971 + randomSeedOffset);
	RandomUniformGenerator lipidProtrusion(-0.035f, 0.035f, 2348971 + randomSeedOffset);
	const float characteristicRadius = std::cbrt(
		midSurfaceRadii.x * midSurfaceRadii.y * midSurfaceRadii.z);
	const std::vector<Float3> surfaceDirections =
		EllipsoidSurfaceDirections(lipidCount, midSurfaceRadii);

	for (const Float3& surfaceDirection : surfaceDirections) {
		const Float3 normal = EllipsoidNormal(surfaceDirection, midSurfaceRadii);
		const Float3 midSurfacePoint = center + surfaceDirection * midSurfaceRadii;
		const float surfaceDisplacement =
			ClosedSurfaceDisplacement(surfaceDirection, characteristicRadius);
		const Float3 lipidCenter = midSurfacePoint + normal
			* (normalOffset + surfaceDisplacement + lipidProtrusion());
		const Float3 lipidDirection = pointOutwards ? normal : -normal;
		const auto [rotationAxis, rotationAngle] = RotationFromZAxisTo(lipidDirection);
		const float selfRotation = randomRotation();
		const Lipids::Select& lipid = getNextRandomLipid();

		std::function<void(Float3&)> transform =
			[lipidCenter, rotationAxis, rotationAngle, selfRotation](Float3& position) {
				position = Float3::rodriguesRotatation(position, Float3{ 0.f, 0.f, 1.f }, selfRotation);
				position = Float3::rodriguesRotatation(position, rotationAxis, rotationAngle);
				position += lipidCenter;
			};
		queuedInsertions.at(lipid.lipidname).emplace_back(
			QueuedInsertion{ *lipid.grofile, std::move(transform), lipid.topfile });
	}
}

static void CreateEllipsoidalMembrane(GroFile& grofile, TopologyFile& topfile,
	const Lipids::Selection& lipidselection, const Float3& center, const Float3& radii,
	const std::string& shapeName, int randomSeed) {
	if (!std::isfinite(center.x) || !std::isfinite(center.y) || !std::isfinite(center.z))
		throw std::invalid_argument(Lima::Format("Membrane {} center must contain finite coordinates.", shapeName));
	if (!std::isfinite(radii.x) || !std::isfinite(radii.y) || !std::isfinite(radii.z)
		|| radii.x <= 0.f || radii.y <= 0.f || radii.z <= 0.f) {
		throw std::invalid_argument(Lima::Format(
			"Membrane {} radii must be positive, finite numbers.", shapeName));
	}

	const float minimumRadius = SimulationBuilder::MinimumSphereRadius(lipidselection);
	const float shortestAxis = std::min({ radii.x, radii.y, radii.z });
	const float longestAxis = std::max({ radii.x, radii.y, radii.z });
	const float minimumCurvatureRadius = shortestAxis * shortestAxis / longestAxis;
	if (minimumCurvatureRadius < minimumRadius) {
		throw std::invalid_argument(Lima::Format(
			"Membrane {} is too tightly curved for the selected lipids; minimum local curvature radius "
			"is {:.3f} nm, but {:.3f} nm is required.",
			shapeName, minimumCurvatureRadius, minimumRadius));
	}

	const float leafletHalfThickness = std::abs(MinParticlePosInDimension(lipidselection, 2));
	const Float3 outerRadii = radii + Float3{ leafletHalfThickness };
	const Float3 innerRadii = radii - Float3{ leafletHalfThickness };
	const int outerLipidCount = std::max(1, static_cast<int>(std::lround(
		EllipsoidSurfaceArea(outerRadii) * lipidDensity)));
	const int innerLipidCount = std::max(1, static_cast<int>(std::lround(
		EllipsoidSurfaceArea(innerRadii) * lipidDensity)));

	std::map<std::string, std::vector<QueuedInsertion>> queuedInsertions;
	for (const auto& lipid : lipidselection)
		queuedInsertions.try_emplace(lipid.lipidname);

	QueueEllipsoidalLeaflet(queuedInsertions, lipidselection, center, radii,
		leafletHalfThickness, outerLipidCount, true, randomSeed);
	QueueEllipsoidalLeaflet(queuedInsertions, lipidselection, center, radii,
		-leafletHalfThickness, innerLipidCount, false, randomSeed + 1);

	int totalIncoming = 0;
	for (const auto& [_, insertions] : queuedInsertions) {
		if (!insertions.empty())
			totalIncoming += static_cast<int>(insertions.front().grofile.atoms.size() * insertions.size());
	}
	grofile.atoms.reserve(grofile.atoms.size() + totalIncoming);

	for (const auto& [_, insertions] : queuedInsertions) {
		for (const QueuedInsertion& insertion : insertions) {
			AddGroAndTopToGroAndTopfile(grofile, insertion.grofile, insertion.positionTransform,
				topfile, insertion.topfile);
		}
	}
}

static void CreateSphericalMembrane(GroFile& grofile, TopologyFile& topfile,
	const Lipids::Selection& lipidselection, const MembraneGeometry::Sphere& sphere, int randomSeed) {
	CreateEllipsoidalMembrane(grofile, topfile, lipidselection, sphere.center,
		Float3{ sphere.radius }, "sphere", randomSeed);
}

static void CreateEllipsoidMembrane(GroFile& grofile, TopologyFile& topfile,
	const Lipids::Selection& lipidselection, const MembraneGeometry::Ellipsoid& ellipsoid, int randomSeed) {
	CreateEllipsoidalMembrane(grofile, topfile, lipidselection, ellipsoid.center,
		ellipsoid.radii, "ellipsoid", randomSeed);
}

void SimulationBuilder::CreateMembrane(GroFile& grofile, TopologyFile& topfile,
	const Lipids::Selection& sharedLipidselection, const MembraneGeometry::Figure& geometry, int randomSeed) {
	validateLipidselection(sharedLipidselection);
	Lipids::Selection lipidselection = sharedLipidselection;
	UseCenteredLipidCopies(lipidselection);
	const std::string name = Lipids::NameSelection(lipidselection);
	if (grofile.title.empty())
		grofile.title = name;
	if (topfile.title.empty())
		topfile.title = name;
	if (!topfile.HasSystem())
		topfile.SetSystem(name);

	std::visit([&](const auto& figure) {
		using FigureType = std::decay_t<decltype(figure)>;
		if constexpr (std::is_same_v<FigureType, MembraneGeometry::Plane>)
			CreatePlanarMembrane(grofile, topfile, lipidselection, figure.z, randomSeed);
		else if constexpr (std::is_same_v<FigureType, MembraneGeometry::Sphere>)
			CreateSphericalMembrane(grofile, topfile, lipidselection, figure, randomSeed);
		else if constexpr (std::is_same_v<FigureType, MembraneGeometry::Ellipsoid>)
			CreateEllipsoidMembrane(grofile, topfile, lipidselection, figure, randomSeed);
	}, geometry);
}

void SimulationBuilder::CreateMembrane(GroFile& grofile, TopologyFile& topfile,
	const Lipids::Selection& lipidselection, float membraneCenter, int randomSeed) {
	CreateMembrane(grofile, topfile, lipidselection, MembraneGeometry::Plane{ membraneCenter }, randomSeed);
}
