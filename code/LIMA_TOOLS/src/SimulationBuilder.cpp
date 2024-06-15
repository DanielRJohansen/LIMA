#include "SimulationBuilder.h"

#include <format>
#include <functional>
#include <algorithm>
#include <random>
#include "BoundaryConditionPublic.h"
#include "EngineCore.h"

const std::array<std::string, 6> LipidSelect::valid_lipids = { "POPC", "POPE", "DDPC", "DMPC", "cholesterol", "DOPC" };


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

Float3 calcDimensions(const GroFile& grofile)
{
	BoundingBox bb{};

	for (const auto& atom : grofile.atoms) {
		for (int dim = 0; dim < 3; dim++) {
			bb.min[dim] = std::min(bb.min[dim], atom.position[dim]);
			bb.max[dim] = std::max(bb.max[dim], atom.position[dim]);
		}		
	}

	const Float3 dims{ bb.max.x - bb.min.x, bb.max.y - bb.min.y, bb.max.z - bb.min.z };
	return dims;
}

float constexpr fursthestDistanceToZAxis(const LipidsSelection& lipidselection) {
	float max_dist = 0;
	for (const auto& lipid : lipidselection) {
		for (const auto& atom : lipid.grofile->atoms) {
			const float dist = sqrtf(atom.position.x * atom.position.x + atom.position.y * atom.position.y);
			max_dist = std::max(max_dist, dist);
		}
	}
	return max_dist;
}



template <typename BondType>
void overwriteBond(const std::vector<BondType>& bonds, std::vector<BondType>& dest, int atomnr_offset) {
	for (const BondType& bond : bonds) {
		dest.emplace_back(bond);
		for (int i = 0; i < bond.n; i++) {
			dest.back().atomGroIds[i] += atomnr_offset;
		}
	}
}

//template <typename BondType>
//void overwriteBond(TopologyFile::SectionRange<BondType> bonds, TopologyFile::SectionRange<BondType> dest, int atomnr_offset) {
//	for (const auto& bond : bonds) {
//		dest.add(bond); // Assuming 'add' appropriately adds items to the container
//		for (int i = 0; i < bond.n; i++) {
//			auto& last = dest.back(); // Reference to modify the last element
//			last.atomGroIds[i] += atomnr_offset;
//		}
//	}
//}


float genRandomAngle() {
	return static_cast<float>(rand() % 360) / 360.f * 2.f * PI;
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

	//outputtopologyfile.atoms.entries.emplace_back(input_atom_top);
	//outputtopologyfile.atoms.entries.back().nr += atom_offset;
	////outputtopologyfile.atoms.entries.back().nr %= 100000;
	//outputtopologyfile.atoms.entries.back().resnr += residue_offset;
	////outputtopologyfile.atoms.entries.back().resnr %= 100000;
}


void constexpr validateLipidselection(const LipidsSelection& lipidselection) {
	int total_percentage = 0;
	for (const auto& lipid : lipidselection) {
		total_percentage += lipid.percentage;
	}
	if (total_percentage != 100) {
		throw std::runtime_error("BuildMembrane failed: Lipid selection did not add up to 100%");
	}

	for (const auto& lipid : lipidselection) {
		if (lipid.grofile->atoms.size() != lipid.topfile->GetLocalAtoms().size()) {
			throw std::runtime_error("BuildMembrane failed: Structure and topology file did not have the same amount of atoms. Please validate your files.");
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

using GetNextRandomLipid = SampleSelectionRandomly<LipidsSelection>;
using GetNextRandomParticle = SampleSelectionRandomly<AtomsSelection>;

namespace SimulationBuilder{

	FilePair buildMembrane(const LipidsSelection& lipidselection, Float3 box_dims) {

		validateLipidselection(lipidselection);
		
		for (auto& lipid : lipidselection) {
			centerMoleculeAroundOrigo(*lipid.grofile);
		}

		const float lipid_density = 1.f / 0.59f;                        // [lipids/nm^2] - Referring to Fig. 6, for DMPC in excess water at 30°C, we find an average cross-sectional area per lipid of A = 59.5 Å2 | https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4241443/
		const float padding = 0.1f;	// [nm]
		//const Float3 mol_dims = calcDimensions(inputgrofile);  // [nm]
		const float molecule_diameter = fursthestDistanceToZAxis(lipidselection) * 2.f;	// [nm]

		const int n_lipids_total = lipid_density * box_dims.x * box_dims.y;
		const Float3 lipids_per_dim_f = sqrtf(static_cast<float>(n_lipids_total));        // We dont xare about z
		Int3 lipids_per_dim{
		static_cast<int>(std::ceil(lipids_per_dim_f.x)),
		static_cast<int>(std::ceil(lipids_per_dim_f.y)),
		1 };
		lipids_per_dim.x += (lipids_per_dim.x % 2) == 0;
		lipids_per_dim.y += (lipids_per_dim.y % 2) == 0;


		const Float3 startpoint = Float3{ box_dims.x * 0.5f };
		const float dist = molecule_diameter + padding;


		auto outputgrofile = std::make_unique<GroFile>();
		outputgrofile->box_size = box_dims;
		outputgrofile->title = "Membrane consisting of ";
		for (const auto& lipid : lipidselection) {
			outputgrofile->title += lipid.lipidname + " (" + std::to_string(lipid.percentage) + "%)    ";
		}
		auto outputtopologyfile = std::make_unique<TopologyFile>();
		outputtopologyfile->name = "monolayer";

		srand(1238971);

		int current_residue_nr_offset = 0;
		GetNextRandomLipid getNextRandomLipid{ lipidselection };

		for (int x = -lipids_per_dim.x / 2; x <= lipids_per_dim.x / 2; x++) {
			for (int y = -lipids_per_dim.y / 2; y <= lipids_per_dim.y / 2; y++) {
				const Float3 center_offset = startpoint + Float3{ static_cast<float>(x), static_cast<float>(y), 0.f } * dist;
				
				const float random_rot = genRandomAngle();
				const int current_atom_nr_offset = outputgrofile->atoms.size();

				const LipidSelect& inputlipid = getNextRandomLipid();
				const auto& inputgrofile = *inputlipid.grofile;
				const auto& inputtopologyfile = *inputlipid.topfile;

				for (int relative_atom_nr = 0; relative_atom_nr < inputgrofile.atoms.size(); relative_atom_nr++) {

					std::function<void(Float3&)> position_transform = [&](Float3& pos) {
						pos.rotateAroundOrigo({ 0.f, 0.f, random_rot });
						pos += center_offset;
						};

					addAtomToFile(*outputgrofile, inputgrofile.atoms[relative_atom_nr], current_atom_nr_offset, current_residue_nr_offset, position_transform);
				}

				outputtopologyfile->AppendTopology(inputlipid.topfile);


				//overwriteBond(inputtopologyfile.singlebonds.entries, outputtopologyfile->singlebonds.entries, current_atom_nr_offset);

				//overwriteBond(inputtopologyfile.GetSinglebonds(), outputtopologyfile->GetSinglebonds(), current_atom_nr_offset);
				//overwriteBond(inputtopologyfile.pairs.entries, outputtopologyfile->pairs.entries, current_atom_nr_offset);
				//overwriteBond(inputtopologyfile.anglebonds.entries, outputtopologyfile->anglebonds.entries, current_atom_nr_offset);
				//overwriteBond(inputtopologyfile.dihedralbonds.entries, outputtopologyfile->dihedralbonds.entries, current_atom_nr_offset);
				//overwriteBond(inputtopologyfile.improperdihedralbonds.entries, outputtopologyfile->improperdihedralbonds.entries, current_atom_nr_offset);

				current_residue_nr_offset++;
			}
		}

		
		return { std::move(outputgrofile), std::move(outputtopologyfile) };
	}

	FilePair makeBilayerFromMonolayer(const FilePair& inputfiles, Float3 box_dims)
	{
		//auto [inputgrofile, inputtopologyfile] = inputfiles;

		const float lowest_zpos = std::min_element(inputfiles.grofile->atoms.begin(), inputfiles.grofile->atoms.end(),
			[](const GroRecord& a, const GroRecord& b) {return a.position.z < b.position.z;}
		)->position.z;	//[nm]

		if (lowest_zpos < 0.f) { throw std::runtime_error("Didnt expect bilayer to go below box"); }


		const float padding_between_layers = 0.05;	// [nm]



		// Copy the existing gro and top file into the output files
		auto outputgrofile = std::make_unique<GroFile>( *inputfiles.grofile );
		outputgrofile->box_size = box_dims;
		outputgrofile->title = "Lipid bi-layer consisting of " + inputfiles.grofile->title;
		auto outputtopologyfile = std::make_unique<TopologyFile>( inputfiles.topfile->copy());
		outputtopologyfile->title = "Lipid bi-layer consisting of " + inputfiles.topfile->title;


		// Translate all the existing postions, so the tails are where the middle of the membrane should be
		const float translation_z = (box_dims.z / 2.f + padding_between_layers/2.f) - lowest_zpos;
		for (auto& atom : outputgrofile->atoms) {
			atom.position.z += translation_z;
		}

		const float lowest_zpos2 = std::min_element(outputgrofile->atoms.begin(), outputgrofile->atoms.end(),
			[](const GroRecord& a, const GroRecord& b) {return a.position.z < b.position.z; }
		)->position.z;	//[nm]

		const int atomnr_offset = inputfiles.grofile->atoms.size();
		const int resnr_offset = inputfiles.grofile->atoms.back().residue_number;

		std::function<void(Float3&)> position_transform = [&](Float3& pos) {
			pos.z = -pos.z;	// Mirror atom in xy plane
			pos.z += box_dims.z - padding_between_layers;	// times 2 since
			//pos.z += lowest_zpos * 2.f - padding;	// Move atom back up to the first layer
			if (pos.z > box_dims.z / 2.f)
				printf("Pos z %f\n", pos.z);
			};

		for (int atom_nr = 0; atom_nr < inputfiles.grofile->atoms.size(); atom_nr++) {
			addAtomToFile(*outputgrofile, outputgrofile->atoms[atom_nr], atomnr_offset, resnr_offset, position_transform);
		}

		for (auto& mol : inputfiles.topfile->molecules.entries) {
			outputtopologyfile->AppendTopology(mol.includeTopologyFile);
		}
		//outputtopologyfile->AppendTopology(inputfiles.topfile);

		//overwriteBond(inputfiles.topfile->GetSinglebonds(), outputtopologyfile->GetSinglebonds(), atomnr_offset);
		//overwriteBond(inputfiles.topfile->singlebonds.entries, outputtopologyfile->singlebonds.entries, atomnr_offset);
		//overwriteBond(inputfiles.topfile->pairs.entries, outputtopologyfile->pairs.entries, atomnr_offset);
		//overwriteBond(inputfiles.topfile->anglebonds.entries, outputtopologyfile->anglebonds.entries, atomnr_offset);
		//overwriteBond(inputfiles.topfile->dihedralbonds.entries, outputtopologyfile->dihedralbonds.entries, atomnr_offset);
		//overwriteBond(inputfiles.topfile->improperdihedralbonds.entries, outputtopologyfile->improperdihedralbonds.entries, atomnr_offset);


		return { std::move(outputgrofile), std::move(outputtopologyfile) };
	}
}

void SimulationBuilder::DistributeParticlesInBox(GroFile& grofile, TopologyFile& topfile, const AtomsSelection& particles, float minDistBetweenAnyParticle, float particlesPerNm3) 
{
	const int blocksPerDim = static_cast<int>(std::ceilf(grofile.box_size.x / 2.f));
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
						topfile.GetLocalAtoms().emplace_back(atomtypeselect.atomtype);
						topfile.GetLocalAtoms().back().id = groId;
						topfile.GetLocalAtoms().back().resnr = resNr;


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