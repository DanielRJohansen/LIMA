#include "SimulationBuilder.h"

#include <format>





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

float fursthestDistanceToZAxis(const std::vector<GroRecord>& atoms) {
	float max_dist = 0;
	for (const auto& atom : atoms) {
		const float dist = sqrtf(atom.position.x * atom.position.x + atom.position.y * atom.position.y);
		max_dist = std::max(max_dist, dist);
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

namespace SimulationBuilder{

	Filepair buildMembrane(Filepair inputfiles, Float3 box_dims) {
		ParsedGroFile inputgrofile = inputfiles.first;
		ParsedTopologyFile inputtopologyfile = inputfiles.second;

		centerMoleculeAroundOrigo(inputgrofile);

		const float lipid_density = 1.f / 0.59f;                        // [lipids/nm^2] - Referring to Fig. 6, for DMPC in excess water at 30°C, we find an average cross-sectional area per lipid of A = 59.5 Å2 | https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4241443/
		//const float lipid_density = 1.f;
		const float padding = 0.1f;	// [nm]
		const Float3 mol_dims = calcDimensions(inputgrofile);  // [nm]
		const float molecule_diameter = fursthestDistanceToZAxis(inputgrofile.atoms) * 2.f;	// [nm]

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
		outputgrofile.title = "Membrane consisting of " + inputgrofile.title;
		ParsedTopologyFile outputtopologyfile{};
		outputtopologyfile.title = "Membrane consisting of " + inputtopologyfile.title;

		srand(1238971);


		int current_residue_nr_offset = 0;

		for (int x = -lipids_per_dim.x / 2; x <= lipids_per_dim.x / 2; x++) {
			for (int y = -lipids_per_dim.y / 2; y <= lipids_per_dim.y / 2; y++) {
				const Float3 center_offset = startpoint + Float3{ static_cast<float>(x), static_cast<float>(y), 0.f } * dist;
				
				const float random_rot = genRandomAngle();
				//const float random_rot = PI/2.f;
				const int current_atom_nr_offset = current_residue_nr_offset * inputgrofile.n_atoms;

				for (int relative_atom_nr = 0; relative_atom_nr < inputgrofile.n_atoms; relative_atom_nr++) {
					outputgrofile.atoms.emplace_back(inputgrofile.atoms[relative_atom_nr]);
					outputgrofile.atoms.back().gro_id += current_atom_nr_offset;
					outputgrofile.atoms.back().gro_id %= 100000;	// Gro_ids only go to 99999
					outputgrofile.atoms.back().position.rotateAroundOrigo({ 0.f, 0.f, random_rot});
					outputgrofile.atoms.back().position += center_offset;
					outputgrofile.atoms.back().residue_number += current_residue_nr_offset;
					outputgrofile.atoms.back().residue_number %= 100000; // Residue ids only go to 99999
					outputgrofile.n_atoms++;

					outputtopologyfile.atoms.entries.emplace_back(inputtopologyfile.atoms.entries[relative_atom_nr]);
					outputtopologyfile.atoms.entries.back().nr += current_atom_nr_offset;
					outputtopologyfile.atoms.entries.back().nr %= 100000;
					outputtopologyfile.atoms.entries.back().resnr += current_residue_nr_offset;
					outputtopologyfile.atoms.entries.back().resnr %= 100000;
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







	
}
