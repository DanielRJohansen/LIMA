// Superstructure to the Filehandler. 
// Provides special functionality for MD files such as .gro .itp .top .pdb and more
#pragma once

#include <filesystem>
#include "Bodies.cuh"
#include "Filehandling.h"

#include <optional>

const bool ENABLE_FILE_CACHING = true;

struct GroRecord {
	int residue_number{};
	std::string residue_name{};
	std::string atom_name{};
	int gro_id{};
	Float3 position{};
	std::optional<Float3> velocity{};
};

struct ParsedGroFile {
	ParsedGroFile() {};
	ParsedGroFile(const std::filesystem::path& path);

	// Contents inside file
	std::string title="";
	std::vector<GroRecord>atoms;
	Float3 box_size{};
	// --------------------- // 

	// Meta data not in file
	std::filesystem::path m_path;
	std::filesystem::file_time_type lastModificationTimestamp;
	bool readFromCache = false;
	// --------------------- // 


	void printToFile() const { printToFile(m_path);};
	void printToFile(const std::filesystem::path& path) const;

	//void addEntry(std::string residue_name, std::string atom_name, const Float3& position);
};

template <typename EntryType>
struct Section {
	std::string title;	// or "directive" in gromacs terms
	const std::string legend;
	std::vector<EntryType> entries;

	std::string composeString() const {
		std::ostringstream oss;

		oss << title << '\n';

		oss << legend << "\n";
		for (const EntryType& entry : entries) {
			entry.composeString(oss);
			oss << '\n';
		}

		oss << '\n';
		return oss.str();
	}
};

class ParsedTopologyFile {	//.top or .itp
public:
	struct MoleculeEntry {
		MoleculeEntry(){}
		MoleculeEntry(const MoleculeEntry& other) : name(other.name) {
			assert(other.includeTopologyFile != nullptr);
			includeTopologyFile = std::make_unique<ParsedTopologyFile>(*other.includeTopologyFile);			
		}
		MoleculeEntry(const std::string& name, std::unique_ptr<ParsedTopologyFile> topfile) 
			: name(name), includeTopologyFile(std::move(topfile)) {}

		std::string name{};	// Name of file without extension and path

		void composeString(std::ostringstream& oss) const;

		// TODO: We need to be able to recursively write these to files too
		std::unique_ptr<ParsedTopologyFile> includeTopologyFile = nullptr;
	};

	
	struct MoleculetypeEntry {
		std::string name;
		int nrexcl{};

		void composeString(std::ostringstream& oss) const;
	};

	// Variable names are from .itp file
	struct AtomsEntry {
		std::optional<std::string> section_name{};// Either a residue or lipid_section

		int nr{};	// Not guaranteed to be unique, atleast not with multiple files!
		std::string type{};
		int resnr{};
		std::string residue{};
		std::string atomname{};
		int cgnr{};
		float charge{};// In elementary charges [e]. Convert to C/mol before using
		float mass{};
		//int chain_id{ -1 };

		void composeString(std::ostringstream& oss) const;
	};

	template <size_t N>
	struct GenericBond {
		static const int n = N;
		int atomGroIds[N]{};
		int funct{};

		void composeString(std::ostringstream& oss) const {
			const int width = 10;
			for (size_t i = 0; i < N; ++i) {
				oss << std::setw(width) << std::right << atomGroIds[i];
			}
			oss << std::setw(width) << std::right << funct;
		}

		void IncrementIds(const int increment)  {
			for (size_t i = 0; i < N; ++i) {
				atomGroIds[i] += increment;
			}
		}	
	};

	using SingleBond = GenericBond<2>;
	using Pair = GenericBond<2>;
	using AngleBond = GenericBond<3>;
	using DihedralBond = GenericBond<4>;
	using ImproperDihedralBond = GenericBond<4>;


	// Create an empty file
	ParsedTopologyFile() {};	
	// Load a file from path
	ParsedTopologyFile(const std::filesystem::path& path);

	void printToFile(const std::filesystem::path& path) const;
	void printToFile() const { printToFile(m_path); };

	// Apply a mapping of resNr and GroID to all entries in this file
	void IncrementIds(int atomNrIncrement, int resNrIncrement);
	ParsedTopologyFile copy() const { return *this; }

	// ----------------------- Information kept in the actual files ----------------------- //
	std::string title="";
	Section<MoleculeEntry> molecules{ "[ molecules ]", generateLegend({}) };
	Section<MoleculetypeEntry> moleculetypes{ "[ moleculetype ]", generateLegend({ "Name", "nrexcl" }) };

	Section<AtomsEntry> atoms{ "[ atoms ]", generateLegend({ "nr", "type", "resnr", "residue", "atom", "cgnr", "charge", "mass" }) };
	Section<SingleBond> singlebonds{ "[ bonds ]", generateLegend({ "ai","aj", "funct","c0","c1","c2","c3" }) };
	Section<Pair> pairs{ "[ pairs ]", generateLegend({ "ai", "aj", "funct", "c0", "c1", "c2", "c3" }) };
	Section<AngleBond> anglebonds{ "[ angles ]", generateLegend({ "ai", "aj", "ak", "funct", "c0", "c1", "c2", "c3" }) };
	Section<DihedralBond> dihedralbonds{ "[ dihedrals ]", generateLegend({ "ai", "aj", "ak", "al", "funct", "c0", "c1", "c2", "c3", "c4", "c5" }) };
	Section<ImproperDihedralBond> improperdihedralbonds{ "[ dihedrals ]", generateLegend({ "ai", "aj", "ak", "al", "funct", "c0", "c1", "c2", "c3" }) };
	// ------------------------------------------------------------------------------------ //


	// ----------------------- Meta data not kept in the file ----------------------- //
	std::filesystem::path m_path;
	std::filesystem::file_time_type lastModificationTimestamp{};
	bool readFromCache = false;
	// ------------------------------------------------------------------------------ //


private:
	// Private copy constructor, since this should be avoided when possible
	ParsedTopologyFile& operator=(const ParsedTopologyFile&) = default;

	static std::string generateLegend(std::vector<std::string> elements);
};



namespace MDFiles {
	namespace fs = std::filesystem;

	struct FilePair {
		std::unique_ptr<ParsedGroFile> grofile;
		std::unique_ptr<ParsedTopologyFile> topfile;
	};

	// Takes the gro and top of "right" and inserts it into left
	void MergeFiles(ParsedGroFile& leftGro, ParsedTopologyFile& leftTop, 
		ParsedGroFile& rightGro, std::unique_ptr<ParsedTopologyFile> rightTop);


	//std::unique_ptr<ParsedTopologyFile> loadTopologyFile(const fs::path& path);



	//struct TrrFile {
	//	static void dumpToFile(const Simulation* sim, const std::string& path);
	//};



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


	// Maybe add simparams here too?
	struct SimulationFilesCollection {
		SimulationFilesCollection() {};
		SimulationFilesCollection(const fs::path& workDir);
		std::unique_ptr<ParsedGroFile> grofile;
		std::unique_ptr<ParsedTopologyFile> topfile;
	};

}