// Superstructure to the Filehandler. 
// Provides special functionality for MD files such as .gro .itp .top .pdb and more
#pragma once

#include <filesystem>
#include "Bodies.cuh"
#include "Filehandling.h"
#include <optional>

struct GroRecord {
	int residue_number{};
	std::string residue_name{};
	std::string atom_name{};
	int gro_id{};
	Float3 position{};
	std::optional<Float3> velocity{};
};

struct ParsedGroFile {
	std::string title;
	int n_atoms{ -1 };
	std::vector<GroRecord>atoms;
	Float3 box_size{};

	void printToFile(const std::filesystem::path& path);
};

struct ParsedTopologyFile {	//.top or .itp
	static std::string generateLegend(std::vector<std::string> elements);

	static const int width = 10;	// Width of elements in files

	struct MoleculeEntry {
		std::string name;
	};
	struct MoleculetypeEntry {
		std::string name;
		int nrexcl{};

		std::string composeString() const;
	};

	// Variable names are from .itp file
	struct AtomsEntry {
		int nr{};	// Not guaranteed to be unique, atleast not with multiple files!
		std::string type{};
		int resnr{};
		std::string residue{};
		std::string atom{};
		int cgnr{};
		float charge{};
		float mass{};
		
		std::string composeString() const;
	};

	template <size_t N>
	struct GenericBond {
		static const int n = N;
		int atom_indexes[N]{};
		int funct{};

		std::string composeString() const {
			std::ostringstream oss;

			for (size_t i = 0; i < N; ++i) {
				oss << std::setw(width) << std::right << atom_indexes[i];
			}
			oss << std::setw(width) << std::right << funct;

			return oss.str();
		}
	};

	using SingleBond = GenericBond<2>;
	using Pair = GenericBond<2>;
	using AngleBond = GenericBond<3>;
	using DihedralBond = GenericBond<4>;
	using ImproperDihedralBond = GenericBond<4>;

	template <typename EntryType>
	struct Section {
		const std::string title;	// or "directive" in gromacs terms
		const std::string legend;
		std::vector<EntryType> entries;

		std::string composeString() {
			std::ostringstream oss;

			oss << title << '\n';

			//oss << EntryType::legend << "\n";
			oss << legend << "\n";
			for (const EntryType& entry : entries) {
				const std::string s = entry.composeString();
				oss << s << '\n';
			}

			oss << '\n';
			return oss.str();
		}
	};


	void printToFile(const std::filesystem::path& path);

	std::string title;

	Section<MoleculeEntry> molecules{ "[ molecule ]", generateLegend({}) };
	Section<MoleculetypeEntry> moleculetypes{ "[ moleculetype ]", generateLegend({ "Name", "nrexcl" }) };

	Section<AtomsEntry> atoms{ "[ atoms ]", generateLegend({ "nr", "type", "resnr", "residue", "atom", "cgnr", "charge", "mass" }) };

	Section<SingleBond> singlebonds{ "[ bonds ]", generateLegend({ "ai","aj", "funct","c0","c1","c2","c3" }) };
	Section<Pair> pairs{ "[ pairs ]", generateLegend({ "ai", "aj", "funct", "c0", "c1", "c2", "c3" }) };
	Section<AngleBond> anglebonds{ "[ angles ]", generateLegend({ "ai", "aj", "ak", "funct", "c0", "c1", "c2", "c3" }) };
	Section<DihedralBond> dihedralbonds{ "[ dihedrals ]", generateLegend({ "ai", "aj", "ak", "al", "funct", "c0", "c1", "c2", "c3", "c4", "c5" }) };
	Section<ImproperDihedralBond> improperdihedralbonds{ "[ dihedrals ]", generateLegend({ "ai", "aj", "ak", "al", "funct", "c0", "c1", "c2", "c3" }) };

	Float3 box_size{};
};


namespace MDFiles {
	namespace fs = std::filesystem;

	ParsedGroFile loadGroFile(const fs::path& path);
	ParsedTopologyFile loadTopologyFile(const fs::path& path);
}