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

float genRandomAngle() {
	return static_cast<float>(rand() % 360) / 360.f * 2.f * PI - PI;
}

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

	outputTopologyFile.AppendMoleculetype(inputTopology->GetMoleculeTypePtr(), inputTopology->forcefieldInclude);
}


void validateLipidselection(const Lipids::Selection& lipidselection) {
	double total_percentage = 0;
	for (const auto& lipid : lipidselection) {
		total_percentage += lipid.percentage;
	}
	if (std::abs(total_percentage - 100) > 0.00001f) {
		throw std::runtime_error("Invalid lipid selection, did not add up to 100%");
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
	SampleSelectionRandomly(const Selection& selection, int seedOffset=0) : selection(selection) {
		// Make a vector that points to the  selection index, so we can randomly select lipids
		for (int i = 0; i < selection.size(); i++) {
			for (int j = 0; j < selection[i].percentage; j++)
				selection_indexes.push_back(i);
		}
		// Shuffle the vector
		g = std::mt19937(34896495 + seedOffset);

		std::shuffle(selection_indexes.begin(), selection_indexes.end(), g);
	}

	const auto& operator()() {
		if (select == 100) {
			select = 0;
			std::shuffle(selection_indexes.begin(), selection_indexes.end(), g);
		}
		return selection[selection_indexes[select++]];
	}

private:
	const Selection& selection;
	std::vector<int> selection_indexes;

	int select = 0;
	std::mt19937 g;
};

using GetNextRandomLipid = SampleSelectionRandomly<Lipids::Selection>;
using GetNextRandomParticle = SampleSelectionRandomly<AtomsSelection>;

void SimulationBuilder::DistributeParticlesInBox(GroFile& grofile, TopologyFile& topfile, const AtomsSelection& particles, float minDistBetweenAnyParticle, float particlesPerNm3) 
{
	const float desiredBlockLen = 2.f;
	const int blocksPerDim = grofile.box_size.x >= 4.f
		? static_cast<int>(std::ceil(grofile.box_size.x / desiredBlockLen))
		: 1;
	const float blockLen = grofile.box_size.x / static_cast<float>(blocksPerDim);

	const float usableBlocklen = blockLen - minDistBetweenAnyParticle;

	const int particlesPerBlock = static_cast<int>(std::ceil(particlesPerNm3 * blockLen * blockLen * blockLen));

	srand(1238971);

	GetNextRandomParticle getNextRandomParticle{ particles };

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

					const Float3 position = Float3{
						static_cast<float>(rand()) / static_cast<float>(RAND_MAX) * usableBlocklen + minDistBetweenAnyParticle * 0.5f,
						static_cast<float>(rand()) / static_cast<float>(RAND_MAX) * usableBlocklen + minDistBetweenAnyParticle * 0.5f,
						static_cast<float>(rand()) / static_cast<float>(RAND_MAX) * usableBlocklen + minDistBetweenAnyParticle * 0.5f
					} + blockStart;

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
						grofile.atoms.emplace_back(GroRecord{ resNr, "XXX", atomtypeselect.atomtype.atomname, groId, position, std::nullopt });

						//grofile.addEntry("XXX", atomtypeselect.atomtype.atomname, position);
						//grofile.atoms.emplace(GroRecord{})

						// Add first the basic atomtype, and then correct the IDs after
						topfile.AppendMolecule(atomtypeselect.atomtype.atomname);
						/*topfile.GetMoleculeType().atoms.emplace_back(atomtypeselect.atomtype);
						topfile.GetMoleculeType().atoms.back().id = groId;
						topfile.GetMoleculeType().atoms.back().resnr = resNr;*/


						//if (topfile.atoms.entries.size() > 1) {
						//	topfile.atoms.entries.back().nr = topfile.atoms.entries[topfile.atoms.entries.size() - 2].nr + 1;
						//	topfile.atoms.entries.back().resnr = topfile.atoms.entries[topfile.atoms.entries.size() - 2].resnr +1;
						//}
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
	int nodesPerDim = 0;

public:
	BoxGrid_(int nodesPerDim) : nodesPerDim(nodesPerDim) {
		nodes.resize(nodesPerDim * nodesPerDim * nodesPerDim);
	}
	
	int GetIndex(NodeIndex nodeindex) {
		return nodeindex.z * nodesPerDim * nodesPerDim + nodeindex.y * nodesPerDim + nodeindex.x;
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
		BoundaryConditionPublic::applyBCNM(absPosHyper, grofile.box_size.x, BoundaryConditionSelect::PBC);

		const NodeIndex nodeindex = NodeIndex{ static_cast<int>(std::floor(absPosHyper.x)), static_cast<int>(std::floor(absPosHyper.y)), static_cast<int>(std::floor(absPosHyper.z)) };

		// Make the positions relative to the nodeIndex
		const Float3 relPos = absPosHyper - Float3{ static_cast<float>(nodeindex.x), static_cast<float>(nodeindex.y), static_cast<float>(nodeindex.z) };
		
		if (relPos.x < 0.f || relPos.y < 0.f || relPos.z < 0.f || relPos.x >= 1.f || relPos.y >= 1.f || relPos.z >= 1.f) {
			throw std::runtime_error("DistributeGrofileparticlesInGrid failed: Particle outside of box");
		}
		

		boxgrid[nodeindex].emplace_back(ParticlePlaceholder{ relPos, true });
	}
}


void SimulationBuilder::SolvateGrofile(GroFile& grofile) {
	if (grofile.box_size.x != ceil(grofile.box_size.x)) {
		throw std::runtime_error("SolvateGroFile failed: Box size must be integers");
	}

	
	const int nodesPerDim = static_cast<int>(grofile.box_size.x);
	int nAtomsInput = grofile.atoms.size();
	BoxGrid_<ParticlePlaceholder> boxgrid{ nodesPerDim };

	DistributeGrofileparticlesInGrid(boxgrid, grofile);

	// TODO: Josiah, is this a problem that our pressure is not precise? If so, we can remove more solvents untill we reach the correct pressure, 
	// but it will be slightly more complex code
	const int desiredSolventsPerNm3 = 34;	// 33.4 is the density of water at 300K, but in some nodes we may have less solvents due to collisions, so we aim a bit higher

	// First add excessive solvents to all blocks
	// TODO: Make OMP
	for (int x = 0; x < nodesPerDim; x++) {
		// The x-column decides the seed
		std::mt19937 rng(x);
		std::uniform_real_distribution<float> dist(0.0f, 1.0f);

		for (int y = 0; y < nodesPerDim; y++) {
			for (int z = 0; z < nodesPerDim; z++) {
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
	for (int x = 0; x < nodesPerDim; x++) {
		for (int y = 0; y < nodesPerDim; y++) {
			for (int z = 0; z < nodesPerDim; z++) {
				const NodeIndex nodeindex = NodeIndex{ x, y, z };
				auto& particles = boxgrid[nodeindex];
				
				// First search through all particles in this block
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

				// Now search through all surrounding blocks
				for (int offsetX = -1; offsetX < 2; offsetX++) {
					for (int offsetY = -1; offsetY < 2; offsetY++) {
						for (int offsetZ = -1; offsetZ < 2; offsetZ++) {
							if (offsetX == 0 && offsetY == 0 && offsetZ == 0)
								continue;

							NodeIndex otherNodeIndex = NodeIndex{ x + offsetX, y + offsetY, z + offsetZ };
							BoundaryConditionPublic::applyBC(otherNodeIndex, nodesPerDim);

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
	
	int countTotal = 0;
	for (int x = 0; x < nodesPerDim; x++) {
		for (int y = 0; y < nodesPerDim; y++) {
			for (int z = 0; z < nodesPerDim; z++) {
				const NodeIndex nodeindex = NodeIndex{ x, y, z };
				int countBlock = 0;
				for (const auto& particle : boxgrid[nodeindex]) {
					if (particle.markedForDeletion || particle.presentInInputfile)
						continue;


					const Float3 absPos = particle.relPos + Float3{ static_cast<float>(x), static_cast<float>(y), static_cast<float>(z) };
					grofile.atoms.push_back(GroRecord{ 1, "SOL", "O", countTotal % 100000, absPos, std::nullopt});

					countTotal++;
					countBlock++;

					if (countBlock == desiredSolventsPerNm3)
						break;
				}
			}
		}
	}
}

void SimulationBuilder::InsertSubmoleculeInSimulation(GroFile& targetGrofile, TopologyFile& targetTopol,
	const GroFile& submolGro, const std::shared_ptr<TopologyFile>& submolTop, Float3 targetCenter)
{
	const Float3 molCenter = MoleculeUtils::GeometricCenter(submolGro);

	std::function<void(Float3&)> position_transform = [&](Float3& pos) {
		pos -= molCenter;
		pos += targetCenter;
	};

	AddGroAndTopToGroAndTopfile(targetGrofile, submolGro, position_transform, targetTopol, submolTop);
}

void SimulationBuilder::InsertSubmoleculesInSimulation(GroFile& targetGrofile, TopologyFile& targetTopol,
	const GroFile& submolGro, const std::shared_ptr<TopologyFile>& submolTop, int nMoleculesToInsert, bool rotateRandomly) 
{
	std::srand(123123123);

	const Float3 molCenter = MoleculeUtils::GeometricCenter(submolGro);
	const float molRadius = MoleculeUtils::Radius(submolGro, molCenter) * 1.1f;

	targetGrofile.atoms.reserve(targetGrofile.atoms.size() + nMoleculesToInsert * submolGro.atoms.size());

	for (int i = 0; i < nMoleculesToInsert; i++) {
		Float3 randomTranslation = Float3{
			static_cast<float>(rand()) / static_cast<float>(RAND_MAX) * (targetGrofile.box_size.x - molRadius*2.f) + molRadius,
			static_cast<float>(rand()) / static_cast<float>(RAND_MAX) * (targetGrofile.box_size.y - molRadius*2.f) + molRadius,
			static_cast<float>(rand()) / static_cast<float>(RAND_MAX) * (targetGrofile.box_size.z - molRadius*2.f) + molRadius
		};
		Float3 randomRotation = Float3{genRandomAngle(), genRandomAngle(), genRandomAngle()};

		std::function<void(Float3&)> position_transform = [&](Float3& pos) {
			pos -= molCenter;
			if (rotateRandomly) {
				Float3::rodriguesRotatation(pos, Float3(0,0,1), randomRotation.z);
				Float3::rodriguesRotatation(pos, Float3(0,1,0), randomRotation.y);
				Float3::rodriguesRotatation(pos, Float3(1,0,0), randomRotation.x);
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
	std::srand(123123123);

	for (auto& lipid : lipidselection) {
		centerMoleculeAroundOrigo(*lipid.grofile);
	}


	GetNextRandomLipid genNextRandomLipid{ lipidselection };


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

		Float3 translationNoise = Float3{
			static_cast<float>(rand()) / static_cast<float>(RAND_MAX) * 0.5,
			static_cast<float>(rand()) / static_cast<float>(RAND_MAX) * 0.5,
			static_cast<float>(rand()) / static_cast<float>(RAND_MAX) * 0.5
		};



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








MDFiles::FilePair SimulationBuilder::CreateMembrane(const Lipids::Selection& lipidselection, Float3 boxSize, float membraneCenter) {
	auto outputgrofile = std::make_unique<GroFile>();
	outputgrofile->box_size = boxSize;
	outputgrofile->title = "Membrane consisting of ";
	for (const auto& lipid : lipidselection) {
		outputgrofile->title += lipid.lipidname + " (" + std::to_string(lipid.percentage) + "%)    ";
	}
	auto outputtopologyfile = std::make_unique<TopologyFile>();
	outputtopologyfile->SetSystem("Membrane");

	CreateMembrane(*outputgrofile, *outputtopologyfile, lipidselection, membraneCenter);

	return { std::move(outputgrofile), std::move(outputtopologyfile) };
}


struct QueuedInsertion {
	GroFile& grofile;
	std::function<void(Float3&)> positionTransform;
	std::shared_ptr<TopologyFile> topfile;
};

void SimulationBuilder::CreateMembrane(GroFile& grofile, TopologyFile& topfile, const Lipids::Selection& lipidselection, float membraneCenter) {

	validateLipidselection(lipidselection);

	for (auto& lipid : lipidselection) {
		centerMoleculeAroundOrigo(*lipid.grofile);
	}

	const float lipid_density = 1.f / 0.59f;                        // [lipids/nm^2] - Referring to Fig. 6, for DMPC in excess water at 30°C, we find an average cross-sectional area per lipid of A = 59.5 Å2 | https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4241443/
	const float lowestZpos = MinParticlePosInDimension(lipidselection, 2);
	const float n_lipids_total = lipid_density * grofile.box_size.x * grofile.box_size.y;
	const int lipidsPerDimx = static_cast<int>(std::ceil(sqrtf(n_lipids_total)));
	const int lipidsPerDimy = static_cast<int>(std::ceil(n_lipids_total / static_cast<float>(lipidsPerDimx)));

	const float distPerX = grofile.box_size.x / static_cast<float>(lipidsPerDimx);
	const float distPerY = grofile.box_size.y / static_cast<float>(lipidsPerDimy);


	const float interLipidLayerSpaceHalf = 0.01f; // [nm]

	srand(1238971);

	GetNextRandomLipid getNextRandomLipid{ lipidselection };

	std::unordered_map<std::string, std::vector<QueuedInsertion>> queuedInsertions;
	for (auto& lipid : lipidselection) {
		queuedInsertions[lipid.lipidname] = std::vector<QueuedInsertion>();
	}

	int nLipidsInserted = 0;
	for (int x = 0; x < lipidsPerDimx; x++) {
		const float packingOffset = x % 2 == 0 ? 0.f : distPerX / 2.f;

		for (int y = 0; y < lipidsPerDimy; y++) {

			if (nLipidsInserted == static_cast<int>(n_lipids_total))
				break;


			const Float3 randomTopDownTranslation{ 0.f,0.f,static_cast<float>(rand()) / static_cast<float>(RAND_MAX) * 0.1f - 0.05f };

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

				//AddGroAndTopToGroAndTopfile(grofile, *inputlipid.grofile, position_transform,
				//	topfile, inputlipid.topfile);

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
