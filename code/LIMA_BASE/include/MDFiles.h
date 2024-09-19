// Superstructure to the Filehandler. 
// Provides special functionality for MD files such as .gro .itp .top .pdb and more
#pragma once

#include <filesystem>
#include "Bodies.cuh"
#include "Filehandling.h"

#include <optional>

#include <stack>
#include <unordered_map>
#include <queue>

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
	fs::file_time_type lastModificationTimestamp;
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
	title, molecules, moleculetype, atoms, bonds, pairs, angles, dihedrals, impropers, position_restraints, _system, cmap, no_section,
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
}







class TopologyFile {	//.top or .itp
public:
	struct ForcefieldInclude;
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
	template <typename EntryType> struct Section {
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


	TopologyFile() {}										// Create an empty file	
	TopologyFile(const fs::path& path);	// Load a file from path	

	void printToFile(const fs::path& path) const;
	void printToFile() const { printToFile(path); };
	void printToFile(const std::string& name) const {
		printToFile(fs::path(path.parent_path() / name));
	}
	void AppendTopology(const std::shared_ptr<TopologyFile>& other);



	// ----------------------- Information kept in the actual files ----------------------- //
	std::string title="";
	Section<MoleculeEntry> molecules{ "[ molecules ]", generateLegend({}) };
	std::optional<Moleculetype> moleculetype;


	std::vector<AtomsEntry>& GetLocalAtoms() { return atoms.entries; }
	const std::vector<AtomsEntry>& GetLocalAtoms() const { return atoms.entries; }
	std::vector<SingleBond>& GetLocalSinglebonds() { return singlebonds.entries; }
	const std::vector<SingleBond>& GetLocalSinglebonds() const { return singlebonds.entries; }
	std::vector<Pair>& GetLocalPairs() { return pairs.entries; }
	const std::vector<Pair>& GetLocalPairs() const { return pairs.entries; }
	std::vector<AngleBond>& GetLocalAnglebonds() { return anglebonds.entries; }
	const std::vector<AngleBond>& GetLocalAnglebonds() const { return anglebonds.entries; }
	std::vector<DihedralBond>& GetLocalDihedralbonds() { return dihedralbonds.entries; }
	const std::vector<DihedralBond>& GetLocalDihedralbonds() const { return dihedralbonds.entries; }
	std::vector<ImproperDihedralBond>& GetLocalImproperDihedralbonds() { return improperdihedralbonds.entries; }
	const std::vector<ImproperDihedralBond>& GetLocalImproperDihedralbonds() const { return improperdihedralbonds.entries; }
	const std::vector<MoleculeEntry>& GetLocalMolecules() const { return molecules.entries; }

	std::vector<ForcefieldInclude> forcefieldIncludes;	// Multiple forcefields can apply to a topology file, in such a case the first forcefield with a hit is used
	std::unordered_map<std::string, std::shared_ptr<TopologyFile>> includeTopologies;
	std::vector<std::string> otherIncludes;
	std::vector<fs::path> GetForcefieldPaths() const;

	template <typename T>
	Section<T>& GetSection() {
		if constexpr (std::is_same_v<T, AtomsEntry>) return atoms;
		else if constexpr (std::is_same_v<T, SingleBond>) return singlebonds;
		else if constexpr (std::is_same_v<T, Pair>) return pairs;
		else if constexpr (std::is_same_v<T, AngleBond>) return anglebonds;
		else if constexpr (std::is_same_v<T, DihedralBond>) return dihedralbonds;
		else if constexpr (std::is_same_v<T, ImproperDihedralBond>) return improperdihedralbonds;
		else if constexpr (std::is_same_v<T, MoleculeEntry>) return molecules;
		else static_assert(std::is_same_v<T, void>, "Unknown section type");
	}

	template <typename T>
	const Section<T>& GetSection() const {
		if constexpr (std::is_same_v<T, AtomsEntry>) return atoms;
		else if constexpr (std::is_same_v<T, SingleBond>) return singlebonds;
		else if constexpr (std::is_same_v<T, Pair>) return pairs;
		else if constexpr (std::is_same_v<T, AngleBond>) return anglebonds;
		else if constexpr (std::is_same_v<T, DihedralBond>) return dihedralbonds;
		else if constexpr (std::is_same_v<T, ImproperDihedralBond>) return improperdihedralbonds;
		else if constexpr (std::is_same_v<T, MoleculeEntry>) return molecules;
		else static_assert(std::is_same_v<T, void>, "Unknown section type");
	}
	// ------------------------------------------------------------------------------------ //

	template<typename T>
	class SectionRange;

	SectionRange<AtomsEntry> GetAllAtoms() const;
	SectionRange<SingleBond> GetAllSinglebonds() const;
	SectionRange<Pair> GetAllPairs() const;
	SectionRange<AngleBond> GetAllAnglebonds() const;
	SectionRange<DihedralBond> GetAllDihedralbonds() const;
	SectionRange<ImproperDihedralBond> GetAllImproperDihedralbonds() const;
	SectionRange<MoleculeEntry> GetAllSubMolecules() const;


	// ----------------------- Meta data not kept in the file ----------------------- //
	std::string name;
	fs::path path;
	fs::file_time_type lastModificationTimestamp{};
	bool readFromCache = false;
	// ------------------------------------------------------------------------------ //

	template <typename T>
	size_t GetElementCount() const {
		size_t count = 0;
		count += GetSection<T>().entries.size();
		for (const auto& mol : molecules.entries) {
			count += mol.includeTopologyFile->GetElementCount<T>();
		}
		return count;
	}


private:
	friend class GenericItpFile;

	static const char commentChar = ';';

	Section<AtomsEntry> atoms{ "[ atoms ]", generateLegend({ "nr", "type", "resnr", "residue", "atom", "cgnr", "charge", "mass" }) };
	Section<SingleBond> singlebonds{ "[ bonds ]", generateLegend({ "ai","aj", "funct","c0","c1","c2","c3" }) };
	Section<Pair> pairs{ "[ pairs ]", generateLegend({ "ai", "aj", "funct", "c0", "c1", "c2", "c3" }) };
	Section<AngleBond> anglebonds{ "[ angles ]", generateLegend({ "ai", "aj", "ak", "funct", "c0", "c1", "c2", "c3" }) };
	Section<DihedralBond> dihedralbonds{ "[ dihedrals ]", generateLegend({ "ai", "aj", "ak", "al", "funct", "c0", "c1", "c2", "c3", "c4", "c5" }) };
	Section<ImproperDihedralBond> improperdihedralbonds{ "[ dihedrals ]", generateLegend({ "ai", "aj", "ak", "al", "funct", "c0", "c1", "c2", "c3" }) };

	// Private copy constructor, since this should be avoided when possible
	TopologyFile& operator=(const TopologyFile&) = default;

	static std::string generateLegend(const std::vector<std::string>& elements);
};

struct TopologyFile::MoleculeEntry {
	MoleculeEntry() {}
	MoleculeEntry(const std::string& name, int globalIndexOfFirstParticle = 0)
		: name(name), globalIndexOfFirstParticle(globalIndexOfFirstParticle) {}
	MoleculeEntry(const std::string& name, std::shared_ptr<TopologyFile> includeTopologyFile,
		int globalIndexOfFirstParticle)
		: name(name), includeTopologyFile(includeTopologyFile),
		globalIndexOfFirstParticle(globalIndexOfFirstParticle) {}

	void composeString(std::ostringstream& oss) const;

	int GlobalIndexOfFinalParticle() const {
		if (includeTopologyFile == nullptr)
			throw std::runtime_error("MoleculeEntry::GlobalIndexOfFinalParticle: includeTopologyFile is nullptr");
		if (includeTopologyFile->atoms.entries.empty())
			throw std::runtime_error("MoleculeEntry::GlobalIndexOfFinalParticle: includeTopologyFile has no atoms");
		return globalIndexOfFirstParticle + includeTopologyFile->atoms.entries.size() - 1;
	}
	int globalIndexOfFirstParticle = 0; // Determined for each molecule when loading the file
	std::string name{};	// Name of file without extension and path
	std::shared_ptr<TopologyFile> includeTopologyFile = nullptr;
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

struct TopologyFile::ForcefieldInclude {
	ForcefieldInclude() {}
	ForcefieldInclude(const fs::path& topolPath, const std::string& includeName);

	/// <summary>
	/// Copies the forcefieldfile AND any include files to a target directory
	/// </summary>
	/// <param name="directory"></param>
	void CopyToDirectory(const fs::path& directory) const;

	bool isUserSupplied;
	fs::path path;
	fs::path name; // Same as path for user supplied, path relative to FF dir for LIMA FF
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
			oss << std::setw(width) << std::right << ids[i];
		}
		oss << std::setw(width) << std::right << funct;
	}

	void IncrementIds(const int increment) {
		for (size_t i = 0; i < N; ++i) {
			ids[i] += increment;
		}
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


template <typename T>
class TopologyFile::SectionRange {
public:
	template <typename TopolPointerType, typename ReturnReferenceType>
	class _Iterator {
	public:
		_Iterator() = default;

		explicit _Iterator(TopolPointerType root)
			: currentTopology(root) {
			moveToNextEntry();
		}

		_Iterator& operator++() {
			moveToNextEntry();
			return *this;
		}

		ReturnReferenceType& operator*()  {
			return currentTopology->template GetSection<T>().entries[elementIndex];
		}

		bool operator==(const _Iterator& other) const {
			return currentTopology == other.currentTopology &&
				elementIndex == other.elementIndex;
		}

		bool operator!=(const _Iterator& other) const {
			return !(*this == other);
		}

	private:
		std::stack<std::queue<TopolPointerType>> topologyStack;
		TopolPointerType currentTopology = nullptr;
		int elementIndex = -1;

		void moveToNextEntry() {
			while (true) {

				// If we have no topology file, go to the next
				if (currentTopology == nullptr) {
					if (topologyStack.empty()) {
						*this = _Iterator();
						return;
					}
					if (topologyStack.top().empty()) {
						topologyStack.pop();
						continue;
					}
					currentTopology = topologyStack.top().front();
					topologyStack.top().pop();
					elementIndex = -1;
				}

				// If our current topology has no more elements, add all include files
				if (elementIndex + 1 == currentTopology->template GetSection<T>().entries.size()) {

					if (currentTopology->molecules.entries.size() > 0) {
						topologyStack.push({});
						for (auto& mol : currentTopology->molecules.entries) {
							topologyStack.top().push(mol.includeTopologyFile.get());
						}

						elementIndex = -1;
						currentTopology = topologyStack.top().front();

						topologyStack.top().pop();
						continue;
					}

					currentTopology = nullptr;
					continue;
				}

				elementIndex++;
				return;
			}
		}
	};
	using Iterator = _Iterator<TopologyFile*, T>;
	using ConstIterator = _Iterator<const TopologyFile*, const T>;

	explicit SectionRange(TopologyFile* const topology) : topology(topology) {
		//topology->LoadAllSubMolecules();
	}

	explicit SectionRange(const TopologyFile* const topology) : topology(const_cast<TopologyFile*>(topology)) {
		//topology->LoadAllTopologies();
	}

	Iterator begin() { return Iterator(topology); }
	Iterator end() { return Iterator(); }
	ConstIterator begin() const { return ConstIterator(topology); }
	ConstIterator end() const { return ConstIterator(); }

	bool operator==(const SectionRange& other) const {
		auto thisIter = begin();
		auto otherIter = other.begin();
		while (thisIter != end() && otherIter != other.end()) {
			const T& thisElement = *thisIter;
			const T& otherElement = *otherIter;

			if (thisElement != otherElement) {
				return false;
			}
			++thisIter;
			++otherIter;
		}
		if (thisIter != end() || otherIter != other.end()) {
			return false;
		}
		return true;
	}


private:
	TopologyFile* const topology = nullptr;
};