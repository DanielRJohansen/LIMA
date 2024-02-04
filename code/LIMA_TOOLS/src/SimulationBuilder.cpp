#include "SimulationBuilder.h"

#include <format>
#include <functional>
#include <algorithm>
#include <random>

const std::array<std::string, 6> LipidSelect::valid_lipids = { "POPC", "POPE", "DDPC", "DMPC", "cholesterol", "DOPC" };


void centerMoleculeAroundOrigo(ParsedGroFile& grofile) {
	if (grofile.n_atoms != grofile.atoms.size()) {
		throw std::runtime_error(std::format("Mismatch between grofiles n_atoms ({}) and actual number of atoms ({})", grofile.n_atoms, grofile.atoms.size()));
	}

	Float3 position_sum{};
	for (const auto& atom : grofile.atoms) {
		position_sum += atom.position;
	}

	const Float3 offset = position_sum / static_cast<float>(grofile.n_atoms);
	for (auto& atom : grofile.atoms) {
		atom.position -= offset;
	}
}

Float3 calcDimensions(const ParsedGroFile& grofile)
{
	BoundingBox bb{};

	for (const auto& atom : grofile.atoms) {
		for (int dim = 0; dim < 3; dim++) {
			*bb.min.placeAt(dim) = std::min(bb.min.at(dim), atom.position.at(dim));
			*bb.max.placeAt(dim) = std::max(bb.max.at(dim), atom.position.at(dim));
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
			dest.back().atom_indexes[i] += atomnr_offset;
		}
	}
}



float genRandomAngle() {
	return static_cast<float>(rand() % 360) / 360.f * 2.f * PI;
}

void addAtomToFile(ParsedGroFile& outputgrofile, ParsedTopologyFile& outputtopologyfile, const GroRecord& input_atom_gro, const ParsedTopologyFile::AtomsEntry input_atom_top, int atom_offset, int residue_offset, 
	std::function<void(Float3&)> position_transform) 
{
	outputgrofile.atoms.emplace_back(input_atom_gro);
	outputgrofile.atoms.back().gro_id += atom_offset;
	outputgrofile.atoms.back().gro_id %= 100000;	// Gro_ids only go to 99999
	outputgrofile.atoms.back().residue_number += residue_offset;
	outputgrofile.atoms.back().residue_number %= 100000; // Residue ids only go to 99999
	position_transform(outputgrofile.atoms.back().position);

	outputgrofile.n_atoms++;

	outputtopologyfile.atoms.entries.emplace_back(input_atom_top);
	outputtopologyfile.atoms.entries.back().nr += atom_offset;
	//outputtopologyfile.atoms.entries.back().nr %= 100000;
	outputtopologyfile.atoms.entries.back().resnr += residue_offset;
	//outputtopologyfile.atoms.entries.back().resnr %= 100000;
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
		if (lipid.grofile->n_atoms != lipid.topfile->atoms.entries.size()) {
			throw std::runtime_error("BuildMembrane failed: Structure and topology file did not have the same amount of atoms. Please validate your files.");
		}
	}
}

struct GetNextRandomLipid {
	GetNextRandomLipid(const LipidsSelection& lipidselection) : lipidselection(lipidselection) {
		// Make a vector that points to the lipid selection index, so we can randomly select lipids
		for (int i = 0; i < lipidselection.size(); i++) {
			for (int j = 0; j < lipidselection[i].percentage; j++)
				lipid_selection_indexes.push_back(i);
		}
		// Shuffle the vector
		g = std::mt19937(34896495);

		std::shuffle(lipid_selection_indexes.begin(), lipid_selection_indexes.end(), g);
	}

	const LipidSelect& operator()() {
		if (lipidselect == 100) {
			lipidselect = 0;
			std::shuffle(lipid_selection_indexes.begin(), lipid_selection_indexes.end(), g);
		}
		return lipidselection[lipid_selection_indexes[lipidselect++]];
	}

private:
	const LipidsSelection& lipidselection;
	std::vector<int> lipid_selection_indexes;

	int lipidselect = 0;
	std::mt19937 g;
};


namespace SimulationBuilder{

	Filepair buildMembrane(const LipidsSelection& lipidselection, Float3 box_dims) {
		//auto [inputgrofile, inputtopologyfile] = inputfiles;

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


		ParsedGroFile outputgrofile{};
		outputgrofile.box_size = box_dims;
		outputgrofile.title = "Membrane consisting of ";
		for (const auto& lipid : lipidselection) {
			outputgrofile.title += lipid.lipidname + " (" + std::to_string(lipid.percentage) + "%)    ";
		}
		ParsedTopologyFile outputtopologyfile{};

		srand(1238971);

		int current_residue_nr_offset = 0;
		GetNextRandomLipid getNextRandomLipid{ lipidselection };

		for (int x = -lipids_per_dim.x / 2; x <= lipids_per_dim.x / 2; x++) {
			for (int y = -lipids_per_dim.y / 2; y <= lipids_per_dim.y / 2; y++) {
				const Float3 center_offset = startpoint + Float3{ static_cast<float>(x), static_cast<float>(y), 0.f } * dist;
				
				const float random_rot = genRandomAngle();
				const int current_atom_nr_offset = outputgrofile.n_atoms;

				const LipidSelect& inputlipid = getNextRandomLipid();
				const auto& inputgrofile = *inputlipid.grofile;
				const auto& inputtopologyfile = *inputlipid.topfile;

				for (int relative_atom_nr = 0; relative_atom_nr < inputgrofile.n_atoms; relative_atom_nr++) {

					std::function<void(Float3&)> position_transform = [&](Float3& pos) {
						pos.rotateAroundOrigo({ 0.f, 0.f, random_rot });
						pos += center_offset;
						};

					addAtomToFile(outputgrofile, outputtopologyfile, inputgrofile.atoms[relative_atom_nr], inputtopologyfile.atoms.entries[relative_atom_nr],
						current_atom_nr_offset, current_residue_nr_offset, position_transform);
				}

				overwriteBond(inputtopologyfile.singlebonds.entries, outputtopologyfile.singlebonds.entries, current_atom_nr_offset);
				overwriteBond(inputtopologyfile.pairs.entries, outputtopologyfile.pairs.entries, current_atom_nr_offset);
				overwriteBond(inputtopologyfile.anglebonds.entries, outputtopologyfile.anglebonds.entries, current_atom_nr_offset);
				overwriteBond(inputtopologyfile.dihedralbonds.entries, outputtopologyfile.dihedralbonds.entries, current_atom_nr_offset);
				overwriteBond(inputtopologyfile.improperdihedralbonds.entries, outputtopologyfile.improperdihedralbonds.entries, current_atom_nr_offset);

				current_residue_nr_offset++;
			}
		}

		
		return { outputgrofile, outputtopologyfile };
	}

	Filepair makeBilayerFromMonolayer(const Filepair& inputfiles, Float3 box_dims)
	{
		auto [inputgrofile, inputtopologyfile] = inputfiles;

		const float lowest_zpos = std::min_element(inputgrofile.atoms.begin(), inputgrofile.atoms.end(), 
			[](const GroRecord& a, const GroRecord& b) {return a.position.z < b.position.z;}
		)->position.z;	//[nm]

		if (lowest_zpos < 0.f) { throw std::runtime_error("Didnt expect bilayer to go below box"); }


		const float padding_between_layers = 0.05;	// [nm]



		// Copy the existing gro and top file into the output files
		ParsedGroFile outputgrofile{ inputgrofile };
		outputgrofile.box_size = box_dims;
		outputgrofile.title = "Lipid bi-layer consisting of " + inputgrofile.title;
		ParsedTopologyFile outputtopologyfile{inputtopologyfile};
		outputtopologyfile.title = "Lipid bi-layer consisting of " + inputtopologyfile.title;


		// Translate all the existing postions, so the tails are where the middle of the membrane should be
		const float translation_z = (box_dims.z / 2.f + padding_between_layers/2.f) - lowest_zpos;
		for (auto& atom : outputgrofile.atoms) {
			atom.position.z += translation_z;
		}

		const float lowest_zpos2 = std::min_element(outputgrofile.atoms.begin(), outputgrofile.atoms.end(),
			[](const GroRecord& a, const GroRecord& b) {return a.position.z < b.position.z; }
		)->position.z;	//[nm]

		const int atomnr_offset = inputgrofile.n_atoms;
		const int resnr_offset = inputgrofile.atoms.back().residue_number;

		std::function<void(Float3&)> position_transform = [&](Float3& pos) {
			pos.z = -pos.z;	// Mirror atom in xy plane
			pos.z += box_dims.z - padding_between_layers;	// times 2 since
			//pos.z += lowest_zpos * 2.f - padding;	// Move atom back up to the first layer
			if (pos.z > box_dims.z / 2.f)
				printf("Pos z %f\n", pos.z);
			};

		for (int atom_nr = 0; atom_nr < inputgrofile.n_atoms; atom_nr++) {
			addAtomToFile(outputgrofile, outputtopologyfile, outputgrofile.atoms[atom_nr], outputtopologyfile.atoms.entries[atom_nr],
				atomnr_offset, resnr_offset, position_transform);
		}

		overwriteBond(inputtopologyfile.singlebonds.entries, outputtopologyfile.singlebonds.entries, atomnr_offset);
		overwriteBond(inputtopologyfile.pairs.entries, outputtopologyfile.pairs.entries, atomnr_offset);
		overwriteBond(inputtopologyfile.anglebonds.entries, outputtopologyfile.anglebonds.entries, atomnr_offset);
		overwriteBond(inputtopologyfile.dihedralbonds.entries, outputtopologyfile.dihedralbonds.entries, atomnr_offset);
		overwriteBond(inputtopologyfile.improperdihedralbonds.entries, outputtopologyfile.improperdihedralbonds.entries, atomnr_offset);


		return { outputgrofile, outputtopologyfile };
	}





	
}
