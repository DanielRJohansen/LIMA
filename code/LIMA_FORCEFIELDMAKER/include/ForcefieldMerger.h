#pragma once

#include "LIMA_FORCEFIELDMAKER/include/ForcefieldTypes.h"
#include "LIMA_FORCEFIELDMAKER/include/Filehandling.h"

#include <assert.h>

class ForcefieldMerger
{
public:
	void mergeForcefields(vector<string> file_paths){

		for (string path : file_paths) {
			vector<vector<string>> rows = Reader::readFile(path, true);
			parseNBAtomtypes(rows, ff_nonbonded);
			parseSinglebondTypes(rows, ff_bondtypes);
			parseAnglebondTypes(rows, ff_angletypes);
			parseDihedralbondTypes(rows, ff_dihedraltypes);
			parseImproperdihedraltypes(rows, ff_improperdihedraltypes);
		}
		printf("\n\n");
		printf("%lld NB_Atomtypes found in param files\n", ff_nonbonded.size());
		printf("%lld Bondtypes found in param files\n", ff_bondtypes.size());
		printf("%lld Angletypes found in param files\n", ff_angletypes.size());
		printf("%lld Dihedraltypes found in param files\n", ff_dihedraltypes.size());
		printf("%lld Improperdihedraltypes found in param files\n", ff_improperdihedraltypes.size());
	}


	vector<NB_Atomtype> ff_nonbonded;
	vector<Singlebondtype> ff_bondtypes;
	vector<Anglebondtype> ff_angletypes;
	vector<Dihedralbondtype> ff_dihedraltypes;
	vector<Improperdihedralbondtype> ff_improperdihedraltypes;

private:
	vector<NB_Atomtype> parseNBAtomtypes(const vector<vector<string>>& rows, vector<NB_Atomtype>& forcefield) {
		STATE current_state = INACTIVE;

		std::map<std::string, float> typeToMass;

		int skipcnt = 0;

		for (const vector<string>& row : rows) {
			if (row.size() == 0)
				continue;

			auto new_state = setState(row[0], current_state);
			if (new_state != current_state && new_state == NONBONDED) {	// Shitty fileformat doesn't use comments around NB
				skipcnt = 2;
			}
			current_state = new_state;

			if (current_state == INACTIVE) { continue; }

			if (skipcnt > 0) {
				skipcnt--;
				continue;
			}

			if (row.size() < 4)
				continue;



			if (current_state == ATOMS) {
				std::string nbtype_name = row[2];
				float mass = stof(row[3]);

				assert(typeToMass.count(nbtype_name) == 0);

				typeToMass.insert(std::pair<std::string, float>(nbtype_name, mass));
			}
			else if (current_state == NONBONDED) {
				const std::string nbtype_name = row[0];
				const float epsilon = abs(stof(row[2]) * 4184.f);			// convert kcal to J								// Do you multiply this together??	// Sometime this value is negative
				const float sigma = stof(row[3]) * 2.f / 10.f;			// convert Rmin/2 to sigma	[A] -> [nm]				// TODO: CHECK WITH ALI :)))))
				const float mass = typeToMass.find(nbtype_name)->second;

				forcefield.push_back(NB_Atomtype{ nbtype_name, mass, sigma, epsilon });

			}

		}

		return forcefield;
	}

	vector<Singlebondtype> parseSinglebondTypes(vector<vector<string>> rows, vector<Singlebondtype>& forcefield) {
		STATE current_state = INACTIVE;

		for (vector<string> row : rows) {
			if (row.size() == 0)
				continue;

			current_state = setState(row[0], current_state);

			if (current_state != BONDS) { continue; }

			if (row.size() < 4) {
				continue;
			}
				
			std::array<string, 2> bonded_typenames{ row[0], row[1] };
			const Singlebondtype bondtype(
				bonded_typenames,
				stof(row[3]) * 0.1f,					// convert A to nm
				stof(row[2]) * 4183.f * 100.f			// convert kcal/(mol*A^2) to J/(mol*nm^2)
			);

			if (isDuplicate(bondtype, forcefield))
				continue;

			forcefield.push_back(bondtype);
		}

		return forcefield;
	}
	vector<Anglebondtype> parseAnglebondTypes(vector<vector<string>> rows, vector<Anglebondtype>& forcefield) {
		STATE current_state = INACTIVE;

		for (vector<string> row : rows) {
			if (row.size() == 0)
				continue;

			current_state = setState(row[0], current_state);

			if (current_state != ANGLES) { continue; }

			if (row.size() < 5)
				continue;

			std::array<string, 3> bonded_typenames{ row[0], row[1], row[2] };
			const Anglebondtype angletype(
				bonded_typenames,
				stof(row[4]) * 2 * 3.1415f / 360.f,					// convert degress to rads
				stof(row[3]) * 4183.f							// convert kcal/(mol*rad^2) to J/(mol*rad^2)
			);
			if (isDuplicate(angletype, forcefield))
				continue;

			forcefield.push_back(angletype);
		}

		return forcefield;
	}
	vector<Dihedralbondtype> parseDihedralbondTypes(vector<vector<string>> rows, vector<Dihedralbondtype>& forcefield) {
		STATE current_state = INACTIVE;

		for (vector<string> row : rows) {
			if (row.size() == 0)
				continue;

			current_state = setState(row[0], current_state);

			if (current_state != DIHEDRALS) {
				continue;
			}

			if (row.size() < 7) {
				continue;
			}
				
			std::array<string, 4> bonded_typenames{ row[0], row[1], row[2], row[3]};
			Dihedralbondtype dihedraltype(
				bonded_typenames,
				stof(row[6]) * 2 * 3.1415f / 360.f,					// convert degress to rads
				stof(row[4]) * 4183.f,							// convert kcal/(mol) to J/(mol)			// Shouldn't this be per rad^2?????
				stoi(row[5])
			);

			if (isDuplicate(dihedraltype, forcefield))
				continue;

			forcefield.push_back(dihedraltype);
		}

		return forcefield;
	}
	vector<Improperdihedralbondtype> parseImproperdihedraltypes(vector<vector<string>> rows, vector<Improperdihedralbondtype>& forcefield) {
		STATE current_state = INACTIVE;

		for (vector<string> row : rows) {
			if (row.size() == 0)
				continue;

			current_state = setState(row[0], current_state);

			if (current_state != IMPROPER) {
				continue;
			}

			if (row.size() < 7) {
				continue;
			}

			std::array<string, 4> bonded_typenames{ row[0], row[1], row[2], row[3] };
			float kpsi = stof(row[4]);	// [kcal/mol/rad^2]
			// ignore row[5]
			float psi0 = stof(row[6]);	// degrees
			Improperdihedralbondtype improperdihedraltype(
				bonded_typenames,
				psi0 * 2.f * PI / 360.f,					// convert degress to rads
				kpsi * 4183.f								// convert kcal/(mol) to J/(mol)
			);

			if (isDuplicate(improperdihedraltype, forcefield))
				continue;

			forcefield.push_back(improperdihedraltype);
		}

		return forcefield;
	}




	// -------------------------- HELPER FUNCTIONS -------------------------- //

	enum STATE { INACTIVE, ATOMS, BONDS, ANGLES, DIHEDRALS, NONBONDED, IMPROPER, CMAP };
	static STATE setState(string s, STATE current_state) {
		if (s == "ATOMS")
			return ATOMS;
		if (s == "BONDS")
			return BONDS;
		if (s == "ANGLES")
			return ANGLES;
		if (s == "DIHEDRALS")
			return DIHEDRALS;
		if (s == "NONBONDED")
			return NONBONDED;
		if (s == "IMPROPER")
			return IMPROPER;
		if (s == "CMAP")
			return CMAP;
		if (s == "NBFIX")
			return INACTIVE;

		return current_state;
	}
	bool isDuplicate(string type, vector<NB_Atomtype>* ff_nonbonded) {
		for (NB_Atomtype atom : *ff_nonbonded) {
			if (atom.type == type)
				return true;
		}
		return false;
	}

	template <typename Bondtype>
	bool isDuplicate(const Bondtype& query_bondtype, const vector<Bondtype>& forcefield_bondtypes) {
		for (const Bondtype& ff_bondtype : forcefield_bondtypes) {
			bool match = 1;

			for (int i = 0; i < Bondtype::n_atoms; i++) {
				match &= (query_bondtype.bonded_typenames[i] == ff_bondtype.bonded_typenames[i]);
			}
			
			if (match) { return true; }
		}
		return false;
	}

};

