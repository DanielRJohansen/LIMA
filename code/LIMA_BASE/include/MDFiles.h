// Superstructure to the Filehandler. 
// Provides special functionality for MD files such as .gro .itp .top .pdb and more
#pragma once

#include <filesystem>
#include "Bodies.cuh"
#include "Filehandling.h"
#include "Simulation.cuh"

#include <optional>
#include <functional>

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
	int n_atoms{ 0 };
	std::vector<GroRecord>atoms;
	Float3 box_size{};

	void printToFile(const std::filesystem::path& path) const;
};

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

struct ParsedTopologyFile {	//.top or .itp

	~ParsedTopologyFile() {
		for (auto& mol : molecules.entries) {
			mol.deletePtr();
		}
	}

	static std::string generateLegend(std::vector<std::string> elements);

	static const int width = 10;	// Width of elements in files

	struct MoleculeEntry {
		std::string name{};

		// TODO: We need to be able to recursively write these to files too

		ParsedTopologyFile* include_file = nullptr;	// TODO: Unsafe, fix with unique_ptr, but i dont know how
		//std::unique_ptr<ParsedTopologyFile> include_files;
		void deletePtr() {
			delete include_file;
		}
	};

	
	struct MoleculetypeEntry {
		std::string name;
		int nrexcl{};

		std::string composeString() const;
	};

	// Variable names are from .itp file
	struct AtomsEntry {
		std::optional<std::string> section_name{};// Either a residue or lipid_section

		int nr{};	// Not guaranteed to be unique, atleast not with multiple files!
		std::string type{};
		int resnr{};
		std::string residue{};
		std::string atom{};
		int cgnr{};
		float charge{};
		float mass{};
		//int chain_id{ -1 };

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


	void printToFile(const std::filesystem::path& path);

	std::string title;

	MoleculeEntry me;
	Section<MoleculeEntry> molecules{ "[ molecule ]", generateLegend({}) };
	Section<MoleculetypeEntry> moleculetypes{ "[ moleculetype ]", generateLegend({ "Name", "nrexcl" }) };

	Section<AtomsEntry> atoms{ "[ atoms ]", generateLegend({ "nr", "type", "resnr", "residue", "atom", "cgnr", "charge", "mass" }) };
	Section<SingleBond> singlebonds{ "[ bonds ]", generateLegend({ "ai","aj", "funct","c0","c1","c2","c3" }) };
	Section<Pair> pairs{ "[ pairs ]", generateLegend({ "ai", "aj", "funct", "c0", "c1", "c2", "c3" }) };
	Section<AngleBond> anglebonds{ "[ angles ]", generateLegend({ "ai", "aj", "ak", "funct", "c0", "c1", "c2", "c3" }) };
	Section<DihedralBond> dihedralbonds{ "[ dihedrals ]", generateLegend({ "ai", "aj", "ak", "al", "funct", "c0", "c1", "c2", "c3", "c4", "c5" }) };
	Section<ImproperDihedralBond> improperdihedralbonds{ "[ dihedrals ]", generateLegend({ "ai", "aj", "ak", "al", "funct", "c0", "c1", "c2", "c3" }) };

	//Float3 box_size{};
};


namespace MDFiles {
	namespace fs = std::filesystem;

	ParsedGroFile loadGroFile(const fs::path& path);
	ParsedGroFile loadGroFile(const Box& box, std::function<Float3(NodeIndex&, Coord&)> getAbsPos);

	std::unique_ptr<ParsedTopologyFile> loadTopologyFile(const fs::path& path);



	struct TrrFile {
		static void dumpToFile(const Simulation* sim, const std::string& path);
	};



	struct ParsedLffFile {
		enum LffSection { title, singlebond, anglebond, dihedralbond, improperdihedralbond, no_section };

		ParsedLffFile(const fs::path& path);
		fs::path path;

		struct Singlebond {
			std::array<uint32_t,2> global_ids;
			float b0;
			float kb;
		};
		struct Anglebond {
			std::array<uint32_t, 3> global_ids;
			float theta0;
			float ktheta;
		};
		struct Dihedralbond {
			std::array<uint32_t, 4> global_ids;
			float phi0;
			float kphi;
			float n;
		};
		struct Improperdihedralbond {
			std::array<uint32_t, 4> global_ids;
			float psi0;
			float kpsi;
		};

		Section<Singlebond> singlebonds;
		Section<Anglebond> anglebonds;
		Section<Dihedralbond> dihedralbonds;
		Section<Improperdihedralbond> improperdihedralbonds;
	};
}