// Superstructure to the FileUtils. 
// Provides special functionality for MD files such as .gro .itp .top .pdb and more
#pragma once

#include <filesystem>
#include "Bodies.cuh"
#include "Filehandling.h"
#include <ranges>
#include <optional>

#include <stack>
#include <unordered_map>
#include <queue>

const bool ENABLE_FILE_CACHING = false;


namespace fs = std::filesystem;

struct GroRecord {
	int residue_number{};
	std::string residueName{};
	std::string atomName{};
	int gro_id{};
	Float3 position{};
	std::optional<Float3> velocity{};
	std::string sourceLine;
};

struct GroFile {
	GroFile() {};
	GroFile(const fs::path& path);

	// Contents inside file
	std::string title="";
	std::vector<GroRecord>atoms;
	Float3 box_size{};
	// --------------------- // 

	// Meta data not in file
	fs::path m_path;
	int64_t lastModificationTimestamp;
	bool readFromCache = false;
	// --------------------- // 


	void printToFile() const { printToFile(m_path);};
	void printToFile(const fs::path& path) const;
	void printToFile(const std::string& name) const {
		printToFile(m_path.parent_path() / name);
	}
};

enum TopologySection {
	// Keywords typically found in topologies
	title, molecules, moleculetype, atoms, bonds, pairs, angles, dihedrals, impropers, position_restraints, _system, // system is a protected keyword??
	cmap, no_section,
	// Keywords typically found in forcefields, but also in topologies when a user wishes to overwrite a parameter
	atomtypes, pairtypes, bondtypes, constainttypes, angletypes, dihedraltypes, impropertypes, defaults, cmaptypes,

	// custom one i defined to do some workarounds
	includes
};

// .itp files follow a simple ruleset. Here we can read such a file, that is NOT a topology file
// This class does NOT support recursive includes, so for reading topologies use the one above
// It is also not very fast, since it stores a lot of data as strings
class GenericItpFile {

	// A section is simply a list of non-comment rows from the file
	using Section = std::vector<std::string>;

	std::unordered_map<TopologySection, Section> sections;
	const Section emptySection;  // Default empty section

public:
	GenericItpFile(const fs::path& path);

	const Section& GetSection(TopologySection section) const {
		auto it = sections.find(section);
		if (it == sections.end()) {
			return emptySection;
		}
		return it->second;
	}
	Section& GetSection(TopologySection section) {
		if (!sections.contains(section)) {
			sections.insert({ section, emptySection });
		}
		return sections.at(section);
	}
};


class PDBfile {
	struct ATOM {
		int atomSerialNumber;
		char atomName[4];
		char altLocIndicator;
		char resName[3];
		char chainID;
		int resSeq;
		char iCode;
		Float3 position; // [nm] (but Angstrom in actual file, so beware of conversion)
		float occupancy;
		float tempFactor;
		char segmentIdentifier[4];
		char elementSymbol;
		char charge[2];
	};
public:

	std::vector<ATOM> ATOMS;
	std::vector<ATOM> HETATMS;
	fs::path mPath;


	PDBfile(const fs::path&);
};

class TopologyFile;

namespace MDFiles {
	namespace fs = std::filesystem;

	struct FilePair {
		std::shared_ptr<GroFile> grofile;
		std::shared_ptr<TopologyFile> topfile;
	};

	// Takes the gro and top of "right" and inserts it into left
	void MergeFiles(GroFile& leftGro, TopologyFile& leftTop,
		GroFile& rightGro, std::shared_ptr<TopologyFile> rightTop);

	// Maybe add simparams here too?
	struct SimulationFilesCollection {
		SimulationFilesCollection() {};
		SimulationFilesCollection(const fs::path& workDir);
		std::unique_ptr<GroFile> grofile;
		std::unique_ptr<TopologyFile> topfile;
	};
	
	struct TrrFile {
		//TrrFile(const fs::path& path);
		TrrFile(Float3 boxSize) : boxSize(boxSize) {};
		void Dump(const fs::path& path) const;
		std::vector<std::vector<Float3>> positions;

	private:
		Float3 boxSize;
	};
}




class TopologyFile {	//.top or .itp
public:
	struct MoleculeEntry;
	struct AtomsEntry;
	template <size_t N>	struct GenericBond;
	struct SingleBond;
	struct Pair;
	struct AngleBond;
	struct DihedralBond;
	struct ImproperDihedralBond;
	struct Moleculetype {
		std::string name{};
		int nrexcl{};
		std::string composeString() const;
	};
	struct Moleculetype1 {
		Moleculetype1() = default;
		Moleculetype1(const TopologyFile& top) {
			*this = top.GetMoleculeType();
		}
		std::string name{};
		int nrexcl{}; // How many consecutive bonds before LJ is enabled again

		std::vector<AtomsEntry> atoms;
		std::vector<SingleBond> singlebonds;
		std::vector<Pair> pairs;
		std::vector<AngleBond> anglebonds;
		std::vector<DihedralBond> dihedralbonds;
		std::vector<ImproperDihedralBond> improperdihedralbonds;		 


		template <typename T>
		std::vector<T>& GetElements() {
			if constexpr (std::is_same_v<T, AtomsEntry>) return atoms;
			else if constexpr (std::is_same_v<T, SingleBond>) return singlebonds;
			else if constexpr (std::is_same_v<T, Pair>) return pairs;
			else if constexpr (std::is_same_v<T, AngleBond>) return anglebonds;
			else if constexpr (std::is_same_v<T, DihedralBond>) return dihedralbonds;
			else if constexpr (std::is_same_v<T, ImproperDihedralBond>) return improperdihedralbonds;
			else static_assert(std::is_same_v<T, void>, "Unknown section type");
		}

		template <typename T>
		const std::vector<T>& GetElements() const {
			if constexpr (std::is_same_v<T, AtomsEntry>) return atoms;
			else if constexpr (std::is_same_v<T, SingleBond>) return singlebonds;
			else if constexpr (std::is_same_v<T, Pair>) return pairs;
			else if constexpr (std::is_same_v<T, AngleBond>) return anglebonds;
			else if constexpr (std::is_same_v<T, DihedralBond>) return dihedralbonds;
			else if constexpr (std::is_same_v<T, ImproperDihedralBond>) return improperdihedralbonds;
			else static_assert(std::is_same_v<T, void>, "Unknown section type");
		}
	};



	struct ForcefieldInclude {
		ForcefieldInclude(const std::string& name = "") : name(name) {};
		ForcefieldInclude(const std::string& name, const fs::path& ownerDir);
		void LoadFullPath(const fs::path& ownerDir); // Necessary when we read from bin files

		/// <summary>
		/// Copies the forcefieldfile AND any include files to a target directory
		/// </summary>
		/// <param name="directory"></param>
		void CopyToDirectory(const fs::path& directory, const fs::path& ownerDir) const;
		const fs::path& Path() const;

		fs::path name; // Either name in resources/forcefields, or a path relative to the topologyfile

	private:
		std::optional<fs::path> path = std::nullopt; // NEVER SAVE/READ THIS TO DISK
	};


	TopologyFile();										// Create an empty file	
	TopologyFile(const fs::path& path, TopologyFile* parentTop=nullptr);	// Load a file from path	

	/// <summary>
	/// Recursively write topology + all includes + forcefieldfiles to the parentpath of the dir
	/// </summary>
	/// <param name="path">path of the output .top file. All other files we be in the same dir</param>
	/// <param name="printForcefieldinclude">Variable used internally in the class</param>
	void printToFile(const fs::path& path, bool printForcefieldinclude=true) const;
	void printToFile() const { printToFile(path); };
	void printToFile(const std::string& name) const {
		printToFile(fs::path(path.parent_path() / name));
	}
	void AppendTopology(const std::shared_ptr<TopologyFile>& other);



	// ----------------------- Information kept in the actual files ----------------------- //
	std::string title="";
	//Section<MoleculeEntry> molecules{ "[ molecules ]", generateLegend({}) };
	//std::optional<Moleculetype> moleculetype;
	struct MoleculeEntry1 {
		std::string name;
		std::shared_ptr<Moleculetype1> moleculetype;
	};

	struct System {
		std::string title;
		std::vector<MoleculeEntry1> molecules;
	};

	const System& GetSystem() const { return system.value();}
	System& GetSystem() { return system.value(); }



	std::optional<ForcefieldInclude> forcefieldInclude = std::nullopt;
	std::unordered_map<std::string, std::shared_ptr<TopologyFile>> includeTopologies;
	std::vector<std::string> otherIncludes;

	template <typename T>
	auto GetAllElements() const {
		// 1. First, transform each MoleculeEntry1 to get the vector of the desired element type.
		// 2. The lambda function does the transformation by calling GetElements<T> on each molecule's Moleculetype1.
		return system->molecules
			| std::views::transform([](const MoleculeEntry1& entry) -> std::vector<T>&{
			return entry.moleculetype->GetElements<T>(); // Retrieve the vector for the specific bond type
				})
			// 3. Join all the individual vectors into a single flattened range of elements (e.g., all SingleBonds from all molecules).
					| std::views::join;
	}


	// Returns the top's forcefield IF it exists.
	// If not, tries to return it's parents IF, if those exists
	// Otherwise, returns nullopt
	std::optional<fs::path> GetForcefieldPath() const;

	const Moleculetype1& GetMoleculeType() const {
		if (moleculetypes1.size() != 1) {
			throw std::runtime_error("Illegal call to GetMoleculeType");
		}
		return *moleculetypes1.begin()->second;
	}
	Moleculetype1& GetMoleculeType() {
		if (moleculetypes1.size() != 1) {
			throw std::runtime_error("Illegal call to GetMoleculeType");
		}
		return *moleculetypes1.begin()->second;
	}


	// ----------------------- Meta data not kept in the file ----------------------- //
	std::string name;
	fs::path path;
	int64_t lastModificationTimestamp;
	bool readFromCache = false;
	// ------------------------------------------------------------------------------ //

	std::unordered_map<std::string, std::shared_ptr<Moleculetype1>> moleculetypes1;
	std::optional<System> system = std::nullopt;

private:
	friend class GenericItpFile;

	static const char commentChar = ';';

	//Section<AtomsEntry> atoms{ "[ atoms ]", generateLegend({ "nr", "type", "resnr", "residue", "atom", "cgnr", "charge", "mass" }) };
	//Section<SingleBond> singlebonds{ "[ bonds ]", generateLegend({ "ai","aj", "funct","c0","c1","c2","c3" }) };
	//Section<Pair> pairs{ "[ pairs ]", generateLegend({ "ai", "aj", "funct", "c0", "c1", "c2", "c3" }) };
	//Section<AngleBond> anglebonds{ "[ angles ]", generateLegend({ "ai", "aj", "ak", "funct", "c0", "c1", "c2", "c3" }) };
	//Section<DihedralBond> dihedralbonds{ "[ dihedrals ]", generateLegend({ "ai", "aj", "ak", "al", "funct", "c0", "c1", "c2", "c3", "c4", "c5" }) };
	//Section<ImproperDihedralBond> improperdihedralbonds{ "[ dihedrals ]", generateLegend({ "ai", "aj", "ak", "al", "funct", "c0", "c1", "c2", "c3" }) };

	//TopologyFile* parentTopology = nullptr; // USE WITH CARE!

	// Private copy constructor, since this should be avoided when possible
	//TopologyFile& operator=(const TopologyFile&) = default;

	static std::string generateLegend(const std::vector<std::string>& elements);

	
	// If this and all submolecules share the same forcefieldinclude, return it, even if this does not have a forcefieldinclude
	// This is useful because GROMACS only allow a single forcefield to be included, and in cases where it is possible
	// we will provide a topology that wont break in GROMACS
	std::optional<ForcefieldInclude> ThisAndAllSubmoleculesShareTheSameForcefieldinclude() const;

	static void ParseFileIntoTopology(TopologyFile&, const fs::path& filepath);

	// Packs atoms and bond information in the moleculetype ptr
	// Returns the next section in the topologyfile
	static TopologySection ParseMoleculetype(std::ifstream& file, std::shared_ptr<Moleculetype1>);
	
};


// Variable names are from .itp file
struct TopologyFile::AtomsEntry {
	std::optional<std::string> section_name{};// Either a residue or lipid_section

	//int nr{};		// Not guaranteed to be unique, atleast not with multiple files!
	int id = -1;	// 0-indexed ID given by LIMA in the order that the atoms are loaded
	std::string type{};	// Atom type
	int resnr{};
	std::string residue{};
	std::string atomname{};	// TODONOW SHould this be prioritized over type???<?>>?>>?<<?!!!!!!
	int cgnr{};
	float charge{};// In elementary charges [e]. Convert to kilo C/mol before using
	float mass{};
	//int chain_id{ -1 };

	void composeString(std::ostringstream& oss) const;

	bool operator==(const AtomsEntry&) const = default;
};



template <size_t N>
struct TopologyFile::GenericBond{
	virtual ~GenericBond() = default;
	static const int n = N;
	//int atomGroIds[N]{};	// We intentionally discard the Incoming id's and give our own ids
	int ids[N]{};	// 0-indexed ID's given by LIMA in the order that the atoms are loaded
	int funct{};

	std::string sourceLine{};	// used for debugging	TODO: remove

	void composeString(std::ostringstream & oss) const {
		const int width = 10;
		for (size_t i = 0; i < N; ++i) {
			oss << std::setw(width) << std::right << ids[i] + 1; // convert back to 1-indexed
		}
		oss << std::setw(width) << std::right << funct;
	}

	bool operator==(const GenericBond<N>& other) const {
		return std::equal(std::begin(ids), std::end(ids), std::begin(other.ids)) && funct == other.funct;
	}
};
struct TopologyFile::SingleBond : GenericBond<2> {};
struct TopologyFile::Pair : GenericBond<2> {};
struct TopologyFile::AngleBond : GenericBond<3> {};
struct TopologyFile::DihedralBond : GenericBond<4> {};
struct TopologyFile::ImproperDihedralBond : GenericBond<4> {};
