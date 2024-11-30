// Superstructure to the FileUtils. 
// Provides special functionality for MD files such as .gro .itp .top .pdb and more
#pragma once

#include "Bodies.cuh"
#include "Filehandling.h"

#include <optional>
#include <filesystem>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <ranges>

const bool ENABLE_FILE_CACHING = true;


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
	int64_t lastModificationTimestamp=0;
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
	includes, notimplemented
};

// .itp files follow a simple ruleset. Here we can read such a file, that is NOT a topology file
// This class does NOT support recursive includes, so for reading topologies use the one above
// It is also not very fast, since it stores a lot of data as strings
class GenericItpFile {

	// A section is simply a list of non-comment rows from the file
	using Section = std::vector<std::string>;

	std::unordered_map<TopologySection, Section> sections;
	Section emptySection;  // Default empty section

public:
	GenericItpFile() {}
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
	struct AtomsEntry;
	template <size_t N, typename ParametersType>	struct GenericBond;
	struct SingleBond;
	struct Pair;
	struct AngleBond;
	struct DihedralBond;
	struct ImproperDihedralBond;
	struct Moleculetype {
		Moleculetype() = default;
		Moleculetype(const std::string& name, int nrexcl) : name(name), nrexcl(nrexcl) {};
		std::string name{};
		int nrexcl{}; // How many consecutive bonds before LJ is enabled again

		std::vector<AtomsEntry> atoms;
		std::vector<SingleBond> singlebonds;
		std::vector<Pair> pairs;
		std::vector<AngleBond> anglebonds;
		std::vector<DihedralBond> dihedralbonds;
		std::vector<ImproperDihedralBond> improperdihedralbonds;		 

		// Only used during parsing!
		std::string mostRecentAtomsSectionName{};
		std::vector<int> groIdToLimaId;

		void ToFile(const fs::path& dir) const;

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
		//ForcefieldInclude(const std::string& name, const fs::path& path) : name(name), path(path) {};
		ForcefieldInclude(const fs::path& filename) : filename(filename) {};

		void SaveToDir(const fs::path& dir) const;
		void AddEntry(TopologySection section, const std::string& entry);

		const fs::path filename; // Either name in resources/forcefields, or a path relative to the topologyfile
		GenericItpFile contents;
	};
	struct MoleculeEntry {
		std::string name{};
		const std::shared_ptr<const Moleculetype> moleculetype = nullptr;
	};
	struct System {
		std::string title{ "noSystem" };
		std::vector<MoleculeEntry> molecules;

		bool IsInit() const { return title != "noSystem"; };
	};

	TopologyFile();										// Create an empty file	
	TopologyFile(const fs::path& path, TopologyFile* parentTop=nullptr);	// Load a file from path	

	/// <summary>
	/// Recursively write topology + all includes + forcefieldfiles to the parentpath of the dir
	/// </summary>
	/// <param name="path">path of the output .top file. All other files we be in the same dir</param>
	/// <param name="printForcefieldinclude">Variable used internally in the class</param>
	void printToFile(const fs::path& path) const;
	void printToFile() const { printToFile(path); };
	void printToFile(const std::string& name) const {
		printToFile(fs::path(path.parent_path() / name));
	}
	



	// ----------------------- Information kept in the actual files ----------------------- //
	std::string title="";







	std::optional<ForcefieldInclude> forcefieldInclude = std::nullopt;
	std::unordered_map<std::string, std::shared_ptr<TopologyFile>> includeTopologies;
	std::vector<std::string> otherIncludes;


	const System& GetSystem() const { 
		if (!m_system.IsInit())
			throw std::runtime_error("System not initialized");
		return m_system;
	}
	void SetSystem(const std::string& systemName) {
		if (m_system.IsInit())
			throw std::runtime_error("System already initialized");
		m_system = System{ systemName };
	}

	template <typename T>
	auto GetAllElements() const {
		// 1. First, transform each MoleculeEntry to get the vector of the desired element type.
		// 2. The lambda function does the transformation by calling GetElements<T> on each molecule's Moleculetype.
		return m_system.molecules
			| std::views::transform(
				[](const MoleculeEntry& entry) -> const std::vector<T>&{return entry.moleculetype->GetElements<T>(); // Retrieve the vector for the specific bond type
				})
			// 3. Join all the individual vectors into a single flattened range of elements (e.g., all SingleBonds from all molecules).
					| std::views::join;
	}

	const GenericItpFile& GetForcefield() const {
		if (forcefieldInclude.has_value()) {
			return forcefieldInclude->contents;
		}
		throw std::runtime_error("No forcefield include in this topology");
	}

	const Moleculetype& GetMoleculeType() const {
		if (moleculetypes.size() != 1) {
			throw std::runtime_error("Illegal call to GetMoleculeType");
		}
		return *moleculetypes.begin()->second;
	}
	Moleculetype& GetMoleculeType() {
		if (moleculetypes.size() != 1) {
			throw std::runtime_error("Illegal call to GetMoleculeType");
		}
		return *moleculetypes.begin()->second;
	}
	const std::shared_ptr<const Moleculetype> GetMoleculeTypePtr() const {
		if (moleculetypes.size() != 1) {
			throw std::runtime_error("Illegal call to GetMoleculeType");
		}
		return moleculetypes.begin()->second;
	}

	// Append a molecule of which the type is already known by the file
	void AppendMolecule(const std::string& moleculename);
	void AppendMoleculetype(const std::shared_ptr<const Moleculetype> moltype, 
		std::optional<ForcefieldInclude> forcefieldInclude=std::nullopt);
	void AppendMolecule(const MoleculeEntry&);
	void AppendMolecules(const std::vector<MoleculeEntry>&);

	// ----------------------- Meta data not kept in the file ----------------------- //
	fs::path path;
	// ------------------------------------------------------------------------------ //

	std::unordered_map<std::string, std::shared_ptr<Moleculetype>> moleculetypes;


private:
	friend class GenericItpFile;

	static const char commentChar = ';';

	std::unordered_set<std::string> defines; // keywords that are defined using #define in a top or itp file. Not fully implemented yet

	/// <summary> Load a .top or .itp into a TopologyFile </summary>
	/// <param name="name">If this is called on an include file, 
	/// this is the name of that include file in the parent file</param>
	static void ParseFileIntoTopology(TopologyFile&, const fs::path& filepath, 
		std::optional<std::string> includefileName =std::nullopt);

	// Packs atoms and bond information in the moleculetype ptr
	// Returns the next section in the topologyfile
	static void ParseMoleculetypeEntry(TopologySection section, 
		const std::string& entry, std::shared_ptr<Moleculetype> moleculetype);

	System m_system{};
};


// Variable names are from .itp file
struct TopologyFile::AtomsEntry {
	std::optional<std::string> section_name{};// Either a residue or lipid_section

	//int nr{};		// Not guaranteed to be unique, atleast not with multiple files!
	int id = -1;	// 0-indexed ID given by LIMA in the order that the atoms are loaded
	std::string type{};	// Atom type
	int resnr{};
	std::string residue{};
	std::string atomname{};
	int cgnr{};
	float charge{};// In elementary charges [e]. Convert to kilo kC/mol before using
	float mass{}; // [g/mol]
	//int chain_id{ -1 };

	void composeString(std::ostringstream& oss) const;

	bool operator==(const AtomsEntry&) const = default;
};



template <size_t N, typename ParametersType>
struct TopologyFile::GenericBond{
	virtual ~GenericBond() = default;
	static const int n = N;
	//int atomGroIds[N]{};	// We intentionally discard the Incoming id's and give our own ids
	std::array<int,N> ids{};	// 0-indexed ID's given by LIMA in the order that the atoms are loaded
	int funct{};

	std::optional<ParametersType> parameters = std::nullopt;

	std::string sourceLine{};	// used for debugging	TODO: remove

	void composeString(std::ostringstream & oss) const {
		const int width = 10;
		for (size_t i = 0; i < N; ++i) {
			oss << std::setw(width) << std::right << ids[i] + 1; // convert back to 1-indexed
		}
		oss << std::setw(width) << std::right << funct << "\n";
	}

    bool operator==(const GenericBond<N, ParametersType>& other) const {
		return std::equal(std::begin(ids), std::end(ids), std::begin(other.ids)) && funct == other.funct;
	}
};
struct TopologyFile::SingleBond : GenericBond<2, Bondtypes::SingleBond::Parameters> {};
struct TopologyFile::Pair : GenericBond<2, Bondtypes::PairBond::Parameters> {};
struct TopologyFile::AngleBond : GenericBond<3, Bondtypes::AngleUreyBradleyBond::Parameters> {};
struct TopologyFile::DihedralBond : GenericBond<4, Bondtypes::DihedralBond::Parameters> {};
struct TopologyFile::ImproperDihedralBond : GenericBond<4, Bondtypes::ImproperDihedralBond::Parameters> {};
