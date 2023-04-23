#pragma once

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <filesystem>
#include <climits>

#include "LIMA_FORCEFIELDMAKER/include/ForcefieldTypes.h"
#include "LIMA_BASE/include/Utilities.h"




namespace fs = std::filesystem;

using namespace std;

// Dunno if this works for folders too
void assertPath(string path) {
	if (!fs::exists(path)) {
		std::cerr << "Could not find path: " << path << "\n";
		abort();
	}
}


class FileHelpers {
public:
	static string pathJoin(string a, string b) { return a + "/" + b; }
};

class Reader {
public:
	static bool ignoreRow(vector<char> ignores, string line) {
		if (line.length() == 0)
			return true;

		char first_char = ' ';
		for (char c : line) {
			if (c != ' ' && c != '\t') {
				first_char = c;
				break;
			}
		}


		for (char c : ignores) {
			if (first_char == c)
				return true;
		}
		return false;
	}


	static vector<string> parseRowExternal(string line) {
		vector<string> row;
		stringstream ss(line);
		string element;
		string e2;

		while (getline(ss, element, ' ')) {
			stringstream ss2 = stringstream(element);
			while (getline(ss2, e2, '\t')) {
				if (e2.length() > 0)
					row.push_back(e2);
			}
		}
		return row;
	}
	static vector <string> parseRowInternal(string line) {
		vector<string> row;
		stringstream ss(line);
		string element;
		string e2;

		while (getline(ss, element, ' ')) {
			stringstream ss2 = stringstream(element);
			while (getline(ss2, e2, ';')) {
				if (e2.length() > 0)
					row.push_back(e2);
			}
		}
		return row;
	}

	static vector<vector<string>> readFile(const string& path, LimaLogger& logger, vector<char> ignores = {';', '#', '!'}, bool internal_file=false) {

		fstream file;
		file.open(path);
		logger.print(std::format("Reading file {}\n", path));

		vector<vector<string>> rows;
		int row_cnt = 0;
		int ignore_cnt = 0;

		string line;
		while (getline(file, line)) {

			if (ignoreRow(ignores, line)) {
				ignore_cnt++;
				continue;
			}

			//cout << line << endl;
			if (internal_file)
				rows.push_back(parseRowInternal(line));
			else
				rows.push_back(parseRowExternal(line));
			row_cnt++;
		}
		//printf("%d rows read. %d rows ignored\n", row_cnt, ignore_cnt);
		logger.print(std::format("{} rows read. {} rows ignored\n", row_cnt, ignore_cnt));


		return rows;
	}
};

class Printer {
public:


	struct FFOutHelpers {
		string includes = "#pragma once\n";

		static string titleH1(string text) {
			return "// ----------------################ " + text + " ################---------------- //\n\n";
		}
		static string titleH2(string text) {
			return "// ---------------- " + text + " ---------------- //\n";
		}
		static string titleH3(string text) {
			return "// " + text + " //\n";
		}
		static string parserTitle(string text) {
			return "# " + text + '\n';
		}
		static string endBlock() {
			return "\n\n\n\n";
		}
	};



	static void printForcefieldSummary(string path, vector<NB_Atomtype> records_nonbonded, Map* map) {
		ofstream file(path, ofstream::out);
		if (!file.is_open()) {
			cout << "Failed to open file " << path << endl;
			exit(0);
		}

		file << FFOutHelpers().includes;
		file << FFOutHelpers::titleH1("Forcefield Non-bonded");
		file << FFOutHelpers::titleH2("Mappings");
		file << FFOutHelpers::parserTitle("mappings");
		for (int i = 0; i < map->n_mappings; i++) {
			file << map->mappings[i].left << ";" << map->mappings[i].right << endl;
		}
		file << FFOutHelpers::endBlock();

		file << FFOutHelpers::titleH2("Non-bonded parameters");
		file << FFOutHelpers::titleH3("{atom_type \t type_id \t mass [g/mol] \t sigma [nm] \t epsilon [J/mol]}");
		file << FFOutHelpers::parserTitle("ff_nonbonded");
		for (NB_Atomtype record : records_nonbonded) {
			file << record.type << ';' << to_string(record.atnum_local) << ';' << to_string(record.mass) << ';' << to_string(record.sigma) << ';' << to_string(record.epsilon) << endl;
		}
		file << FFOutHelpers::endBlock();


		file.close();
	}

	static void printForcefield(string path, vector<Atom> atoms, vector<Bondtype> bonds, vector<Angletype> angles, vector<Dihedraltype> dihedrals) {

		ofstream file(path, ofstream::out);
		if (!file.is_open()) {
			printf(("Failed to open file\n"));
			exit(0);
		}

		file << FFOutHelpers::titleH1("Forcefield Non-bonded");
		file << FFOutHelpers::titleH3("Atoms {particle id [simulation specific]\tatomtype_id [simulation specific]}");
		file << FFOutHelpers::parserTitle("atoms");
		for (Atom atom : atoms) {
			file << to_string(atom.id) << ";" << to_string(atom.atomtype_id) << endl;
		}
		file << FFOutHelpers::endBlock();






		file << FFOutHelpers::titleH2("Bonds");
		file << FFOutHelpers::titleH3("{ID_p1 \t ID_p2 \t Atomtype \t Atomtype \t b_0 [nm] \t k_b [J/(mol * nm^2)]}");
		file << FFOutHelpers::parserTitle("bonds");
		for (Bondtype bond : bonds) {
			file << to_string(bond.id1) << ';' << to_string(bond.id2) << ';'
				<< bond.type1 << ';' << bond.type2 << ';'
				<< to_string(bond.b0) << ';' << to_string(bond.kb) << endl;
		}
		file << FFOutHelpers::endBlock();



		file << FFOutHelpers::titleH2("Angles");
		file << FFOutHelpers::titleH3("{Atom-IDs \t Atomtypes \t theta_0 [rad] \t k_theta [J/(mol * rad^2)}");
		file << FFOutHelpers::parserTitle("angles");
		for (Angletype angle : angles) {
			file << to_string(angle.id1) << ';' << to_string(angle.id2) << ';' << to_string(angle.id3) << ';'
				<< angle.type1 << ';' << angle.type2 << ';' << angle.type3 << ';'
				<< to_string(angle.theta0) << ';' << to_string(angle.ktheta) << endl;
		}
		file << FFOutHelpers::endBlock();



		file << FFOutHelpers::titleH2("Dihedrals");
		file << FFOutHelpers::titleH3("{Atom IDs \t Atomtypes \t phi_0 [rad] \t k_phi [J/(mol * rad^2)] \t n}");
		file << FFOutHelpers::parserTitle("dihedrals");
		for (Dihedraltype dihedral : dihedrals) {
			file << to_string(dihedral.id1) << ';' << to_string(dihedral.id2) << ';' << to_string(dihedral.id3) << ';' << to_string(dihedral.id4) << ';'
				<< dihedral.type1 << ';' << dihedral.type2 << ';' << dihedral.type3 << ';' << dihedral.type4 << ';'
				<< to_string(dihedral.phi0) << ';' << to_string(dihedral.kphi) << ';' << to_string(dihedral.n) << endl;
		}
		file << FFOutHelpers::endBlock();

		file.close();
	}

	static void printFFNonbonded(string path, vector<NB_Atomtype> forcefield) {
		ofstream file(path, ofstream::out);
		if (!file.is_open()) {
			cout << "Failed to open file " << path << endl;
			exit(0);
		}
		file << FFOutHelpers::titleH1("Forcefield Non-bonded");
		file << FFOutHelpers::titleH2("Non-bonded parameters");
		file << FFOutHelpers::titleH3("{atom_type \t mass [g/mol] \t sigma [nm] \t epsilon [J/mol]}");
		file << FFOutHelpers::parserTitle("ff_nonbonded");
		for (NB_Atomtype atomtype : forcefield) {
			file << atomtype.type << ';' << to_string(atomtype.mass) << ';' << to_string(atomtype.sigma) << ';' << to_string(atomtype.epsilon) << endl;
		}
		file << FFOutHelpers::endBlock();

		file.close();

	}

	static void printFFBonded(string path, vector<Bondtype> bondtypes, vector<Angletype> angletypes, vector<Dihedraltype> dihedraltypes) {
		ofstream file(path, ofstream::out);
		if (!file.is_open()) {
			cout << "Failed to open file " << path << endl;
			exit(0);
		}


		file << FFOutHelpers::titleH1("Forcefield Bonded");


		file << FFOutHelpers::titleH2("Bondtype parameters");
		file << FFOutHelpers::titleH3("{atom_types \t b_0 [nm] \t k_b [J/(mol*nm^2)] \t }");
		file << FFOutHelpers::parserTitle("ff_bondtypes");
		for (Bondtype bondtype : bondtypes) {
			file << bondtype.type1 << ';' << bondtype.type2 << ';' << to_string(bondtype.b0) << ';' << to_string(bondtype.kb) << endl;
		}
		file << FFOutHelpers::endBlock();



		file << FFOutHelpers::titleH2("Angletype parameters");
		file << FFOutHelpers::titleH3("{atom_types \t theta_0 [rad] \t k_theta [J/(mol*rad^2)] \t }");
		file << FFOutHelpers::parserTitle("ff_angletypes");
		for (Angletype angle : angletypes) {
			file << angle.type1 << ';' << angle.type2 << ';' << angle.type3 << ';' << to_string(angle.theta0) << ';' << to_string(angle.ktheta) << endl;
		}
		file << FFOutHelpers::endBlock();



		file << FFOutHelpers::titleH2("Dihedraltype parameters");
		file << FFOutHelpers::titleH3("{atom_types \t phi_0 [rad] \t k_phi [J/(mol)] \t n}");
		file << FFOutHelpers::parserTitle("ff_dihedraltypes");
		for (Dihedraltype dihedral : dihedraltypes) {
			file << dihedral.type1 << ';' << dihedral.type2 << ';' << dihedral.type3 << ';' << dihedral.type4 << ';' << to_string(dihedral.phi0) << ';' << to_string(dihedral.kphi) << ';' << to_string(dihedral.n) << endl;
		}
		file << FFOutHelpers::endBlock();



		file.close();
	}
};



