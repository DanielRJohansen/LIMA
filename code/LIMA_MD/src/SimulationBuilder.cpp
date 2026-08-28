#include "SimulationBuilder.h"

#include "BoundaryConditionPublic.h"
#include "EngineCore.h"
#include "Statistics.h"
#include "MoleculeGraph.h"
#include "MoleculeUtils.h"

#include <format>
#include <functional>
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

float constexpr fursthestDistanceToZAxis(const Lipids::Selection& lipidselection) {
	float max_dist = 0;
	for (const auto& lipid : lipidselection) {
		for (const auto& atom : lipid.grofile->atoms) {
			const float dist = sqrtf(atom.position.x * atom.position.x + atom.position.y * atom.position.y);
			max_dist = std::max(max_dist, dist);
		}
	}
	return max_dist;
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
		throw std::runtime_error(std::format("Invalid lipid selection, did not add up to 100% {:.2f}", total_percentage));
	}

	for (const auto& lipid : lipidselection) {
		//if (lipid.grofile->atoms.size() != lipid.topfile->GetLocalAtoms().size()) {
		if (lipid.grofile->atoms.size() != lipid.topfile->GetMoleculeType().atoms.size()) {
			throw std::runtime_error(std::format("BuildMembrane failed: Structure and topology file did not have the same amount of atoms. Please validate your files.\nGRO:{}\nTOP:{}",
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




struct ParticlePlaceholder {
	Float3 relPos{};	// [nm]
	bool presentInInputfile = false;	// We need to know the difference between the particles already present, and the ones we just added
	bool markedForDeletion = false;
};

// TODO: remove this temp class

template<typename T>
class BoxGrid_ {	// TODO: Rename
	std::vector<std::vector<T>> nodes;
	Int3 nodesPerDim = 0;

public:
	BoxGrid_(Int3 nodesPerDim) : nodesPerDim(nodesPerDim) {
		nodes.resize(nodesPerDim.x * nodesPerDim.y * nodesPerDim.z);
	}
	
	int GetIndex(const NodeIndex& nodeindex) {
		return nodeindex.z * nodesPerDim.x * nodesPerDim.y + nodeindex.y * nodesPerDim.x + nodeindex.x;
	}

	std::vector<T>& operator[](NodeIndex index3d) {
		BoundaryConditionPublic::applyBC(index3d, nodesPerDim);
		return nodes[GetIndex(index3d)];
	}
	const std::vector<T>& operator[](NodeIndex index3d) const {
		BoundaryConditionPublic::applyBC(index3d, nodesPerDim);
		return nodes[GetIndex(index3d)];
	}
};



void DistributeGrofileparticlesInGrid(BoxGrid_<ParticlePlaceholder>& boxgrid, const GroFile& grofile) {
	for (const auto& atom : grofile.atoms) {
		Float3 absPosHyper = atom.position;
		BoundaryConditionPublic::applyBCNM(absPosHyper, grofile.box_size, BoundaryConditionSelect::PBC);

		const NodeIndex nodeindex = NodeIndex{ static_cast<int>(std::floor(absPosHyper.x)), static_cast<int>(std::floor(absPosHyper.y)), static_cast<int>(std::floor(absPosHyper.z)) };

		// Make the positions relative to the nodeIndex
		const Float3 relPos = absPosHyper - Float3{ static_cast<float>(nodeindex.x), static_cast<float>(nodeindex.y), static_cast<float>(nodeindex.z) };
		
		if (relPos.x < 0.f || relPos.y < 0.f || relPos.z < 0.f || relPos.x >= 1.f || relPos.y >= 1.f || relPos.z >= 1.f) {
			throw std::runtime_error("DistributeGrofileparticlesInGrid failed: Particle outside of box");
		}
		

		boxgrid[nodeindex].emplace_back(ParticlePlaceholder{ relPos, true });
	}
}


void SimulationBuilder::SolvateGrofile(GroFile& grofile, TopologyFile& topfile, int desiredSolventsPerNm3) {

	throw std::runtime_error("SolvateGrofile is not implemented yet");
	if (grofile.box_size.x != ceil(grofile.box_size.x)) {
		throw std::runtime_error("SolvateGroFile failed: Box size must be integers");
	}

	
	const Int3 gridDim = grofile.box_size.ToInt3();
	const int nAtomsInput = grofile.atoms.size();
	const int startResidueId = grofile.atoms.empty() ? 1 : grofile.atoms.back().residue_number + 1;

	BoxGrid_<ParticlePlaceholder> boxgrid{ gridDim };


	//BoxGrid_<Float3> nonSolventPositions{ gridDim };	// This is used to store the positions of the particles that are not solvents, so we can remove them later
	//for (auto elem : grofile.atoms) {				
	//	NodeIndex nodeindex = NodeIndex{ static_cast<int>(std::floor(elem.position.x)), static_cast<int>(std::floor(elem.position.y)), static_cast<int>(std::floor(elem.position.z)) };
	//	Float3 relpos = elem.position - Float3{ static_cast<float>(nodeindex.x), static_cast<float>(nodeindex.y), static_cast<float>(nodeindex.z) };

	//	nonSolventPositions[nodeindex].emplace_back(relpos);
	//}

	DistributeGrofileparticlesInGrid(boxgrid, grofile);

	// TODO: Josiah, is this a problem that our pressure is not precise? If so, we can remove more solvents untill we reach the correct pressure, 
	// but it will be slightly more complex code

	// First add excessive solvents to all blocks
	// TODO: Make OMP
	for (int x = 0; x < gridDim.x; x++) {
		// The x-column decides the seed
		std::mt19937 rng(x);
		std::uniform_real_distribution<float> dist(0.0f, 1.0f);

		for (int y = 0; y < gridDim.y; y++) {
			for (int z = 0; z < gridDim.z; z++) {
				const NodeIndex nodeindex = NodeIndex{ x, y, z };
				auto& particles = boxgrid[nodeindex];
				for (int i = 0; i < desiredSolventsPerNm3 + 20; i++) {	// +20 so we can remove any particles that are too close
					const Float3 relPos = Float3{ dist(rng), dist(rng), dist(rng) };
					particles.emplace_back(ParticlePlaceholder{ relPos, false });
				}
			}
		}
	}


	const float distanceThreshold = 0.12;	// [nm]

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


						if ((particles[pid].relPos - particles[otherId].relPos).len() < distanceThreshold) {
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

							NodeIndex otherNodeIndex = NodeIndex{ x + offsetX, y + offsetY, z + offsetZ };
							BoundaryConditionPublic::applyBC(otherNodeIndex, gridDim);

							const Float3 relPosOffsetOther = Float3{ static_cast<float>(offsetX), static_cast<float>(offsetY), static_cast<float>(offsetZ) };

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

									if ((particle.relPos - queryParticlePosition).len() < distanceThreshold) {
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
	
	const float bondLen = 0.1f;	// probably shouldnt be hardcoded here...
	const float bondAngle = 109.47f * PI / 180.f;	// [rad] 
	const Float3 _h1Pos = Float3::rodriguesRotatation(Float3(0.f, 0.f, -bondLen), Float3(0,1,0), -bondAngle * 0.5f);
	const Float3 _h2Pos = Float3::rodriguesRotatation(Float3(0.f, 0.f, -bondLen), Float3(0, 1, 0), bondAngle * 0.5f);

	//float zOffset = -h1Pos.z * 0.5f;
	RandomUniformGeneratorUnitvector genRandomUnitVector(1238971);	// Seed offset so we can get different random vectors for each simulation
	RandomUniformGenerator genRandomAngle(-PI, PI);

	int atomCount = nAtomsInput;
	int solventCount = 0;
	for (int x = 0; x < gridDim.x; x++) {
		for (int y = 0; y < gridDim.y; y++) {
			for (int z = 0; z < gridDim.z; z++) {
				const NodeIndex nodeindex = NodeIndex{ x, y, z };
				int nSolventsInBlock = 0;
				for (const auto& solvent : boxgrid[nodeindex]) {
					if (solvent.markedForDeletion || solvent.presentInInputfile)
						continue;


					Float3 rotVector = genRandomUnitVector();
					float rotAngle = genRandomAngle();

					Float3 h1Pos = Float3::rodriguesRotatation(_h1Pos, rotVector, rotAngle);
					Float3 h2Pos = Float3::rodriguesRotatation(_h2Pos, rotVector, rotAngle);

					const Float3 blockOffset = Float3{ static_cast<float>(x), static_cast<float>(y), static_cast<float>(z) };
					grofile.atoms.push_back(GroRecord{ (solventCount+startResidueId) % 100000, SmallString("SOL"), SmallString("OW"),  (atomCount + 1)% 100000, solvent.relPos + blockOffset, std::nullopt});
					grofile.atoms.push_back(GroRecord{ (solventCount+startResidueId) % 100000, SmallString("SOL"), SmallString("HW1"), (atomCount + 2)% 100000, solvent.relPos + blockOffset + h1Pos, std::nullopt });
					grofile.atoms.push_back(GroRecord{ (solventCount+startResidueId) % 100000, SmallString("SOL"), SmallString("HW2"), (atomCount + 3)% 100000, solvent.relPos + blockOffset + h2Pos, std::nullopt });

					nSolventsInBlock++;
					atomCount += 3;
					solventCount++;
					if (nSolventsInBlock >= desiredSolventsPerNm3)
						break;
				}
			}
		}
	}

	//topfile.AppendSolvents(solventCount, FileUtils::GetLimaDir() / "resources" / "forcefields" / "charmm27.ff" / "spce.itp");
	//topfile.AppendSolvents()
	/*TopologyFile solventTop{ FileUtils::GetLimaDir() / "resources" / "forcefields" / "charmm27.ff" / "spce.itp" };
	topfile.AppendMoleculetype(solventTop.GetMoleculeTypePtr(), solventTop.forcefieldInclude);
	for (size_t i = 0; i < solventCount; i++) {
		topfile.AppendMolecule("SOL");
	}*/
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

	for (auto& lipid : lipidselection) {
		centerMoleculeAroundOrigo(*lipid.grofile);
	}


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








MDFiles::FilePair SimulationBuilder::CreateMembrane(const Lipids::Selection& lipidselection, Float3 boxSize,
	const MembraneGeometry::Figure& geometry) {
	auto outputgrofile = std::make_unique<GroFile>();
	outputgrofile->box_size = boxSize;
	outputgrofile->title = "Membrane consisting of ";
	for (const auto& lipid : lipidselection) {
		outputgrofile->title += lipid.lipidname + " (" + std::to_string(lipid.percentage) + "%)    ";
	}
	auto outputtopologyfile = std::make_unique<TopologyFile>();
	outputtopologyfile->SetSystem("Membrane");

	CreateMembrane(*outputgrofile, *outputtopologyfile, lipidselection, geometry);

	return { std::move(outputgrofile), std::move(outputtopologyfile) };
}

MDFiles::FilePair SimulationBuilder::CreateMembrane(const Lipids::Selection& lipidselection, Float3 boxSize,
	float membraneCenter) {
	return CreateMembrane(lipidselection, boxSize, MembraneGeometry::Plane{ membraneCenter });
}


struct QueuedInsertion {
	GroFile& grofile;
	std::function<void(Float3&)> positionTransform;
	std::shared_ptr<TopologyFile> topfile;
};

static void CreatePlanarMembrane(GroFile& grofile, TopologyFile& topfile,
	const Lipids::Selection& lipidselection, float membraneCenter) {

	const float lowestZpos = MinParticlePosInDimension(lipidselection, 2);
	const float n_lipids_total = lipidDensity * grofile.box_size.x * grofile.box_size.y; // (per side)
	const int lipidsPerDimx = static_cast<int>(std::ceil(sqrtf(n_lipids_total)));
	const int lipidsPerDimy = static_cast<int>(std::ceil(n_lipids_total / static_cast<float>(lipidsPerDimx)));

	const float distPerX = grofile.box_size.x / static_cast<float>(lipidsPerDimx);
	const float distPerY = grofile.box_size.y / static_cast<float>(lipidsPerDimy);


	const float interLipidLayerSpaceHalf = 0.01f; // [nm]

	RandomUniformGenerator genRandomAngle(-PI, PI);
	GetNextRandomLipid getNextRandomLipid{ lipidselection };
	RandomUniformGenerator genRandomUpDownTranslation(-0.05f, 0.05f);

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


			const Float3 randomTopDownTranslation{ 0.f,0.f,genRandomUpDownTranslation()};

			// Insert top lipid
			{
				const Lipids::Select& inputlipid = getNextRandomLipid();

				const Float3 lipidCenter = Float3{
					static_cast<float>(x) * distPerX + distPerX / 2.f,
					static_cast<float>(y) * distPerY + distPerY / 4.f + packingOffset,
					membraneCenter + std::abs(lowestZpos) + interLipidLayerSpaceHalf
				};


				const float randomRot = genRandomAngle();
				std::function<void(Float3&)> position_transform = [lipidCenter, randomRot, randomTopDownTranslation](Float3& pos) {
					pos = Float3::rodriguesRotatation(pos, Float3{ 0,0,1 }, randomRot);
					pos += lipidCenter;
					pos += randomTopDownTranslation;
					};

				queuedInsertions.at(inputlipid.lipidname).emplace_back(QueuedInsertion{ *inputlipid.grofile, position_transform, inputlipid.topfile });
			}

			// Insert bottom lipid
			{
				const Lipids::Select& inputlipid = getNextRandomLipid();

				const Float3 lipidCenter = Float3{
					static_cast<float>(x) * distPerX + distPerX / 2.f,
					static_cast<float>(y) * distPerY + distPerY / 4.f + packingOffset,
					membraneCenter - std::abs(lowestZpos) - interLipidLayerSpaceHalf
				};

				const float randomRot = genRandomAngle();
				std::function<void(Float3&)> position_transform = [lipidCenter, randomRot, randomTopDownTranslation](Float3& pos) {
					pos = Float3::rodriguesRotatation(pos, Float3{ 0,0,1 }, randomRot);
					pos = Float3::rodriguesRotatation(pos, Float3{ 1,0,0 }, PI); // Rotate 180 degrees around x-axis
					pos += lipidCenter;
					pos += randomTopDownTranslation;
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

static void QueueSphericalLeaflet(
	std::map<std::string, std::vector<QueuedInsertion>>& queuedInsertions,
	const Lipids::Selection& lipidselection,
	const Float3& sphereCenter,
	float lipidCenterRadius,
	int lipidCount,
	bool pointOutwards,
	int randomSeedOffset)
{
	GetNextRandomLipid getNextRandomLipid{ lipidselection, randomSeedOffset };
	RandomUniformGenerator randomRotation(-PI, PI, 1238971 + randomSeedOffset);
	RandomUniformGenerator radialNoise(-0.05f, 0.05f, 2348971 + randomSeedOffset);
	const float goldenAngle = std::numbers::pi_v<float> * (3.f - std::sqrt(5.f));

	for (int i = 0; i < lipidCount; ++i) {
		const float z = 1.f - 2.f * (static_cast<float>(i) + 0.5f) / static_cast<float>(lipidCount);
		const float xyRadius = std::sqrt(std::max(0.f, 1.f - z * z));
		const float theta = goldenAngle * static_cast<float>(i);
		const Float3 radialDirection{
			xyRadius * std::cos(theta),
			xyRadius * std::sin(theta),
			z
		};

		const Float3 lipidCenter = sphereCenter
			+ radialDirection * (lipidCenterRadius + radialNoise());
		const Float3 lipidDirection = pointOutwards ? radialDirection : -radialDirection;
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

static void CreateSphericalMembrane(GroFile& grofile, TopologyFile& topfile,
	const Lipids::Selection& lipidselection, const MembraneGeometry::Sphere& sphere) {
	if (!std::isfinite(sphere.radius) || sphere.radius <= 0.f)
		throw std::invalid_argument("Membrane sphere radius must be a positive, finite number.");
	if (!std::isfinite(sphere.center.x) || !std::isfinite(sphere.center.y) || !std::isfinite(sphere.center.z))
		throw std::invalid_argument("Membrane sphere center must contain finite coordinates.");

	const float minimumRadius = SimulationBuilder::MinimumSphereRadius(lipidselection);
	if (sphere.radius < minimumRadius) {
		throw std::invalid_argument(std::format(
			"Membrane sphere radius {:.3f} nm is too small for the selected lipids; minimum is {:.3f} nm.",
			sphere.radius, minimumRadius));
	}

	const float leafletHalfThickness = std::abs(MinParticlePosInDimension(lipidselection, 2));
	const float outerLipidCenterRadius = sphere.radius + leafletHalfThickness;
	const float innerLipidCenterRadius = sphere.radius - leafletHalfThickness;
	const int outerLipidCount = std::max(1, static_cast<int>(std::lround(
		4.f * std::numbers::pi_v<float> * outerLipidCenterRadius * outerLipidCenterRadius * lipidDensity)));
	const int innerLipidCount = std::max(1, static_cast<int>(std::lround(
		4.f * std::numbers::pi_v<float> * innerLipidCenterRadius * innerLipidCenterRadius * lipidDensity)));

	std::map<std::string, std::vector<QueuedInsertion>> queuedInsertions;
	for (const auto& lipid : lipidselection)
		queuedInsertions.try_emplace(lipid.lipidname);

	QueueSphericalLeaflet(queuedInsertions, lipidselection, sphere.center,
		outerLipidCenterRadius, outerLipidCount, true, 0);
	QueueSphericalLeaflet(queuedInsertions, lipidselection, sphere.center,
		innerLipidCenterRadius, innerLipidCount, false, 1);

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

void SimulationBuilder::CreateMembrane(GroFile& grofile, TopologyFile& topfile,
	const Lipids::Selection& lipidselection, const MembraneGeometry::Figure& geometry) {
	validateLipidselection(lipidselection);
	for (const auto& lipid : lipidselection)
		centerMoleculeAroundOrigo(*lipid.grofile);

	std::visit([&](const auto& figure) {
		using FigureType = std::decay_t<decltype(figure)>;
		if constexpr (std::is_same_v<FigureType, MembraneGeometry::Plane>)
			CreatePlanarMembrane(grofile, topfile, lipidselection, figure.z);
		else if constexpr (std::is_same_v<FigureType, MembraneGeometry::Sphere>)
			CreateSphericalMembrane(grofile, topfile, lipidselection, figure);
	}, geometry);
}

void SimulationBuilder::CreateMembrane(GroFile& grofile, TopologyFile& topfile,
	const Lipids::Selection& lipidselection, float membraneCenter) {
	CreateMembrane(grofile, topfile, lipidselection, MembraneGeometry::Plane{ membraneCenter });
}
