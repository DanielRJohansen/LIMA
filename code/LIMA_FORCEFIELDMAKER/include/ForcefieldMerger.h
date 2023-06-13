#pragma once

#include "LIMA_FORCEFIELDMAKER/include/ForcefieldTypes.h"
#include "LIMA_BASE/include/Filehandling.h"


class ForcefieldMerger
{
public:
	void mergeForcefields(vector<string> file_paths, LimaLogger& logger){

		for (string path : file_paths) {
			vector<vector<string>> rows = Reader::readFile(path, logger);
			parseNBAtomtypes(rows, &ff_nonbonded);
			parseBondtypes(rows, &ff_bondtypes);
			parseAngletypes(rows, &ff_angletypes);
			parseDihedraltypes(rows, &ff_dihedraltypes);
		}
		printf("\n\n");
		printf("%lld NB_Atomtypes found in param files\n", ff_nonbonded.size());
		printf("%lld Bondtypes found in param files\n", ff_bondtypes.size());
		printf("%lld Angletypes found in param files\n", ff_angletypes.size());
		printf("%lld Dihedraltypes found in param files\n", ff_dihedraltypes.size());
	}


	vector<NB_Atomtype> ff_nonbonded;
	vector<Singlebondtype> ff_bondtypes;
	vector<Anglebondtype> ff_angletypes;
	vector<Dihedralbondtype> ff_dihedraltypes;

private:

	void parseNBAtomtypes(vector<vector<string>> rows, vector<NB_Atomtype>* ff_nonbonded) {
		STATE current_state = INACTIVE;
		NB_Atomtype* nb_atomtype;


		for (vector<string> row : rows) {
			if (row.size() == 0)
				continue;

			current_state = setState(row[0], current_state);

			switch (current_state) {
			case ATOMS:
				if (row.size() < 4)
					continue;

				if (isDuplicate(row[2], ff_nonbonded))
					continue;

				ff_nonbonded->push_back(NB_Atomtype(row[2], stof(row[3])));
				break;

			case NONBONDED:
				if (row.size() < 4)
					continue;

				nb_atomtype = getPtrToNBAtomtype(row[0], ff_nonbonded);
				if (nb_atomtype == nullptr)						// In this case the atom was not found
					break;

				nb_atomtype->epsilon = abs(stof(row[2]) * 4184.f);			// convert kcal to J								// Do you multiply this together??	// Sometime this value is negative
				nb_atomtype->sigma = stof(row[3]) * 2.f / 10.f;			// convert Rmin/2 to sigma	[A] -> [nm]				// TODO: CHECK WITH ALI :)))))
				break;

			default:
				break;
			}
		}
	}

	void parseBondtypes(vector<vector<string>> rows, vector<Singlebondtype>* ff_bondtypes) {
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
			Singlebondtype bondtype(
				bonded_typenames,
				stof(row[3]) * 0.1f,					// convert A to nm
				stof(row[2]) * 4183.f * 100.f			// convert kcal/(mol*A^2) to J/(mol*nm^2)
			);

			if (isDuplicate(bondtype, ff_bondtypes))
				continue;

			ff_bondtypes->push_back(bondtype);

		}
	}
	void parseAngletypes(vector<vector<string>> rows, vector<Anglebondtype>* ff_angletypes) {
		STATE current_state = INACTIVE;


		for (vector<string> row : rows) {
			if (row.size() == 0)
				continue;

			current_state = setState(row[0], current_state);

			if (current_state != ANGLES) { continue; }


			if (row.size() < 5)
				continue;

			std::array<string, 3> bonded_typenames{ row[0], row[1], row[2] };
			Anglebondtype angletype(
				bonded_typenames,
				stof(row[4]) * 2 * 3.1415f / 360.f,					// convert degress to rads
				stof(row[3]) * 4183.f							// convert kcal/(mol*rad^2) to J/(mol*rad^2)
			);
			if (isDuplicate(angletype, ff_angletypes))
				continue;

			ff_angletypes->push_back(angletype);

		}
	}
	void parseDihedraltypes(vector<vector<string>> rows, vector<Dihedralbondtype>* ff_dihedraltypes) {			// what the fuck is multiplicity?!??!!?!?!?
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

			if (isDuplicate(dihedraltype, ff_dihedraltypes))
				continue;

			ff_dihedraltypes->push_back(dihedraltype);

		}
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
		return current_state;
	}
	bool isDuplicate(string type, vector<NB_Atomtype>* ff_nonbonded) {
		for (NB_Atomtype atom : *ff_nonbonded) {
			if (atom.type == type)
				return true;
		}
		return false;
	}
	bool isDuplicate(Singlebondtype bondtype, vector<Singlebondtype>* ff_bondtypes) {
		for (Singlebondtype bt : *ff_bondtypes) {
			if (bondtype.bonded_typenames[0] == bt.bonded_typenames[0] && bondtype.bonded_typenames[1] == bt.bonded_typenames[1])
				return true;
		}
		return false;
	}
	bool isDuplicate(Anglebondtype angletype, vector<Anglebondtype>* ff_angletypes) {
		for (Anglebondtype at : *ff_angletypes) {
			if (angletype.bonded_typenames[0] == at.bonded_typenames[0] && angletype.bonded_typenames[1] == at.bonded_typenames[1] && angletype.bonded_typenames[2] == at.bonded_typenames[2])
				return true;
		}
		return false;
	}
	bool isDuplicate(Dihedralbondtype dihedraltype, vector<Dihedralbondtype>* ff_dihedraltypes) {
		for (Dihedralbondtype dt : *ff_dihedraltypes) {
			if (dihedraltype.bonded_typenames[0] == dt.bonded_typenames[0] 
				&& dihedraltype.bonded_typenames[1] == dt.bonded_typenames[1] 
				&& dihedraltype.bonded_typenames[2] == dt.bonded_typenames[2] 
				&& dihedraltype.bonded_typenames[3] == dt.bonded_typenames[3])
				return true;
		}
		return false;
	}

	NB_Atomtype* getPtrToNBAtomtype(string type, vector<NB_Atomtype>* ff_nonbonded) {
		for (int i = 0; i < ff_nonbonded->size(); i++) {
			NB_Atomtype* atom = &(ff_nonbonded->at(i));
			if (atom->type == type)
				return atom;
		}
		return nullptr;
	}

};

