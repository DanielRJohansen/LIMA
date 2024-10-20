#include "MDFiles.h"
#include "Filehandling.h"
#include "MDFilesSerialization.h"

#include <algorithm>
#include <format>


using namespace FileUtils;
using namespace MDFiles;
namespace fs = std::filesystem;

GroRecord parseGroLine(const std::string& line) {
	GroRecord record;

	// Parse residue number (5 positions, integer) directly
	record.residue_number = std::stoi(line.substr(0, 5));

	// Direct assignment avoiding unnecessary erase-remove idiom
	// Use std::string::find_first_not_of and find_last_not_of to trim spaces
	auto trimSpaces = [](const std::string& str) -> std::string {
		size_t start = str.find_first_not_of(' ');
		size_t end = str.find_last_not_of(' ');
		return start == std::string::npos ? "" : str.substr(start, end - start + 1);
		};

	// Parse residue name (5 characters)
	record.residueName = trimSpaces(line.substr(5, 5));

	// Parse atom name (5 characters)
	record.atomName = trimSpaces(line.substr(10, 5));

	// Parse atom number (5 positions, integer) directly
	record.gro_id = std::stoi(line.substr(15, 5));

	// Parse position (in nm, x y z in 3 columns, each 8 positions with 3 decimal places)
	// Using std::strtod for direct parsing from substring to double
	auto parseDouble = [](const std::string& str) -> double {
		return std::strtod(str.c_str(), nullptr);
		};

	record.position.x = parseDouble(line.substr(20, 8));
	record.position.y = parseDouble(line.substr(28, 8));
	record.position.z = parseDouble(line.substr(36, 8));

	if (line.size() >= 68) {
		record.velocity = Float3{};
		// Parse velocity (in nm/ps (or km/s), x y z in 3 columns, each 8 positions with 4 decimal places)
		record.velocity->x = parseDouble(line.substr(44, 8));
		record.velocity->y = parseDouble(line.substr(52, 8));
		record.velocity->z = parseDouble(line.substr(60, 8));
	}

	return record;
}

std::string composeGroLine(const GroRecord& record) {
	std::ostringstream oss;

	// Format and write each part of the GroRecord
	oss << std::setw(5) << std::left << record.residue_number
		<< std::setw(5) << std::left << record.residueName
		<< std::setw(5) << std::left << record.atomName
		<< std::setw(5) << std::right << record.gro_id
		<< std::setw(8) << std::fixed << std::setprecision(3) << record.position.x
		<< std::setw(8) << std::fixed << std::setprecision(3) << record.position.y
		<< std::setw(8) << std::fixed << std::setprecision(3) << record.position.z;

	// Append velocity if present
	if (record.velocity) {
		oss << std::setw(8) << std::fixed << std::setprecision(4) << record.velocity->x
			<< std::setw(8) << std::fixed << std::setprecision(4) << record.velocity->y
			<< std::setw(8) << std::fixed << std::setprecision(4) << record.velocity->z;
	}

	return oss.str();
}

//void TopologyFile::MoleculeEntry::composeString(std::ostringstream& oss) const {
//	throw std::runtime_error("Dont call this function");
//	//if (includeTopologyFile->name == "")
//	//	throw std::runtime_error("Trying to save a topology with an empty include topology name");
//	//oss << includeTopologyFile->name;// << " 1";
//}

std::string TopologyFile::Moleculetype::composeString() const {
	std::ostringstream oss;
	const std::string legend = generateLegend({ "Name", "nrexcl" });
	oss << "[ moleculetype ]\n";
	oss << legend << "\n";
	oss << std::right << std::setw(10) << name << std::setw(10) << nrexcl << "\n";
	return oss.str();
}

void TopologyFile::AtomsEntry::composeString(std::ostringstream& oss) const {
	if (section_name) {
		oss << section_name.value() << "\n";
	}
	oss << std::right
		<< std::setw(10) << id + 1 // convert back to 1-indexed
		<< std::setw(10) << type
		<< std::setw(10) << resnr
		<< std::setw(10) << residue
		<< std::setw(10) << atomname
		<< std::setw(10) << cgnr
		<< std::setw(10) << std::fixed << std::setprecision(2) << charge
		<< std::setw(10) << std::fixed << std::setprecision(3) << mass;
}


GroFile::GroFile(const fs::path& path) : m_path(path){
	if (!(path.extension().string() == std::string{ ".gro" }))
		throw std::runtime_error("Expected .gro extension");
	if (!fs::exists(path))
		throw std::runtime_error(std::format("File \"{}\" was not found", path.string()));

	lastModificationTimestamp = TimeSinceEpoch(fs::last_write_time(path));

	if (UseCachedBinaryFile(path)) {
		readGroFileFromBinaryCache(path, *this);
	}
	else {
		assert(path.extension() == ".gro");
		if (!fs::exists(path)) { throw std::runtime_error(std::format("File \"{}\" was not found", path.string())); }

		std::ifstream file;
		file.open(path);
		if (!file.is_open() || file.fail()) {
			throw std::runtime_error(std::format("Failed to open file {}\n", path.string()).c_str());
		}

		int skipCnt = 2;	// First 2 lines are title and atom count

		const int min_chars = 5 + 5 + 5 + 5 + 8 + 8 + 8;
		int nAtoms = 0;
		// Forward declaring for optimization reasons
		std::string line{}, word{}, prevLine{};
		while (getline(file, line)) {

			if (skipCnt > 0) {
				if (skipCnt == 2) {
					// 1st line is title
					title = line;
				}
				if (skipCnt == 1) {
					// 2nd line is atom count
					nAtoms = std::stoi(line);
					atoms.reserve(nAtoms);
				}

				skipCnt--;
				continue;
			}


			if (prevLine.empty()) {
				prevLine = line;
				continue;
			}

			if (line.empty()) {
				assert(false);	// I guess this shouldn't happen?
			}
			else {
				if (!(prevLine.length() >= min_chars))
					int c = 0;
				assert(prevLine.length() >= min_chars);
				atoms.emplace_back(parseGroLine(prevLine));
				atoms.back().sourceLine = prevLine;
			}
			prevLine = line;
		}

		if (!prevLine.empty()) {
			int dim = 0;
			std::stringstream ss(prevLine);
			while (std::getline(ss, word, ' ')) {
				if (!word.empty()) {
					box_size[dim++] = std::stof(word);
					if (dim == 3)
						break;
				}
			}
		}

		assert(atoms.size() == nAtoms);

		// Save a binary cached version of the file to we can read it faster next time
		WriteFileToBinaryCache(*this);
	}
}

//
//class TopologySectionGetter {
//	int dihedralCount = 0;
//	int dihedraltypesCount = 0;
//public:
//	TopologySection operator()(const std::string& directive) {
//		if (directive == "molecules") return molecules;
//		if (directive == "moleculetype") return moleculetype;
//		if (directive == "atoms") return atoms;
//		if (directive == "bonds") return bonds;
//		if (directive == "pairs") return pairs;
//		if (directive == "angles") return angles;
//		if (directive == "dihedrals") {
//			if (++dihedralCount <= 2) return dihedralCount == 1 ? dihedrals : impropers;
//			throw std::runtime_error("Encountered 'dihedral' directive more than 2 times in .itp/.top file");
//		}
//		if (directive == "position_restraints") return position_restraints;
//		if (directive == "system") return _system;
//		if (directive == "cmap") return cmap;
//
//		if (directive == "atomtypes") return atomtypes;
//		if (directive == "pairtypes") return pairtypes;
//		if (directive == "bondtypes") return bondtypes;
//		if (directive == "constrainttypes") return constainttypes;
//		if (directive == "angletypes") return angletypes;
//		if (directive == "dihedraltypes") {
//			if (++dihedraltypesCount <= 2) return dihedraltypesCount == 1 ? dihedraltypes : impropertypes;
//			throw std::runtime_error("Encountered 'dihedraltypes' directive more than 2 times in .itp/ file");
//		}
//		if (directive == "impropertypes") return impropertypes;
//		if (directive == "defaults") return defaults;
//		if (directive == "cmaptypes") return cmaptypes;
//
//		throw std::runtime_error(std::format("Got unexpected topology directive: {}", directive));
//	}
//};
//
//std::string extractSectionName(const std::string& line) {
//	size_t start = line.find('[');
//	size_t end = line.find(']', start);
//	if (start != std::string::npos && end != std::string::npos) {
//		// Extract the text between '[' and ']'
//		std::string sectionName = line.substr(start + 1, end - start - 1);
//		removeWhitespace(sectionName);
//		return sectionName;
//	}
//	return "";
//}
//
//
//std::string GetCleanFilename(const fs::path& path) {
//	auto filename = path.stem().string();
//	const std::string prefix = "topol_";
//	return filename.starts_with(prefix) ? filename.substr(prefix.length()) : filename;
//}
//
//std::optional<fs::path> _SearchForFile(const fs::path& dir, const std::string& filename) {
//	const std::array<std::string, 2> extensions = { ".itp", ".top" };
//	const std::array<std::string, 2> prefixes = { std::string(""), "topol_" };
//
//	for (const auto& ext : extensions) {
//		for (const auto& prefix : prefixes) {
//			fs::path path = dir / (prefix + filename + ext);
//			if (fs::exists(path)) {
//				return path;
//			}
//		}
//	}
//
//	return std::nullopt;
//}
//
//fs::path SearchForFile(const fs::path& workDir, const std::string& includeName) {
//	// First look relative to the current topology file
//	std::optional<fs::path> includePath = _SearchForFile(workDir, includeName);
//
//	// If no value, look in the default includes dir
//	if (!includePath.has_value())
//		includePath = _SearchForFile(FileUtils::GetLimaDir() / "resources/Slipids", includeName);
//
//	if (!includePath.has_value())
//		throw std::runtime_error(std::format("Could not find file \"{}\" in directory \"{}\"", includeName, workDir.string()));
//
//	return includePath.value();
//}
//
//
//
//template <int n>
//bool VerifyAllParticlesInBondExists(const std::vector<int>& groIdToLimaId, int ids[n]) {
//	for (int i = 0; i < n; i++) {
//		if (ids[i] >= groIdToLimaId.size() || groIdToLimaId[ids[i]] == -1)
//			return false;
//	}
//	return true;
//}
//
///// <returns>True if we change section/stay on no section, in which case we should 'continue' to the next line. False otherwise</returns>
//bool HandleTopologySectionStartAndStop(const std::string& line, TopologySection& currentSection, TopologySectionGetter& sectionGetter) {
//
//	if (line.empty() && currentSection == TopologySection::title) {
//		currentSection = no_section;
//		return true;
//	}
//	else if (!line.empty() && line[0] == '[') {
//		currentSection = sectionGetter(extractSectionName(line));
//		return true;
//	}
//	
//	return false; // We did not change section, and should continue parsing the current line
//}
//bool isOnlySpacesAndTabs(const std::string& str) {
//	return std::all_of(str.begin(), str.end(), [](char c) {
//		return c == ' ' || c == '\t';
//		});
//}

void GroFile::printToFile(const std::filesystem::path& path) const {
	if (path.extension().string() != ".gro") { throw std::runtime_error(std::format("Got {} extension, expected .gro", path.extension().string())); }

	std::ofstream file(path);
	if (!file.is_open()) {
		throw std::runtime_error(std::format("Failed to open file {}", path.string()));
	}

	// Print the title and number of atoms
	file << title << "\n";
	file << atoms.size() << "\n";

	// Iterate over atoms and print them
	for (const auto& atom : atoms) {
		// You need to define how GroRecord is formatted
		file << composeGroLine(atom) << "\n";
	}

	// Print the box size
	file << box_size.x << " " << box_size.y << " " << box_size.z << "\n";

	file.close();

	// Also cache the file
	WriteFileToBinaryCache(*this, path);
}








std::string TopologyFile::generateLegend(const std::vector<std::string>& elements)
{
	std::ostringstream legend;
	legend << ';'; // Start with a semicolon

	for (const auto& element : elements) {
		legend << std::setw(10) << std::right << element;
	}
	return legend.str();
}



void TopologyFile::printToFile(const std::filesystem::path& path, bool printForcefieldinclude) const {
	//const auto ext = path.extension().string();
	//if (ext != ".top" && ext != ".itp") { throw std::runtime_error(std::format("Got {} extension, expectec [.top/.itp]", ext)); }

	//{
	//	std::ofstream file(path);
	//	if (!file.is_open()) {
	//		throw std::runtime_error(std::format("Failed to open file {}", path.string()));
	//	}
	//	
	//	file << title << "\n\n";

	//	if (printForcefieldinclude) {
	//		std::optional<ForcefieldInclude> sharedff = ThisAndAllSubmoleculesShareTheSameForcefieldinclude();
	//		if (sharedff.has_value()) {
	//			sharedff.value().CopyToDirectory(path.parent_path(), path.parent_path());
	//			file << ("#include \"" + sharedff->name.string() + "\"\n");
	//			file << "\n";
	//			printForcefieldinclude = false;
	//		}
	//		else if (forcefieldInclude) {
	//			forcefieldInclude.value().CopyToDirectory(path.parent_path(), path.parent_path());
	//			file << ("#include \"" + forcefieldInclude->name.string() + "\"\n");
	//			file << "\n";
	//		}
	//	}

	//	for (const auto& include : includeTopologies) {
	//		file << ("#include \"topol_" + include.first + ".itp\"\n");
	//	}
	//	file << "\n";

	//	if (moleculetype.has_value()) { file << moleculetype.value().composeString(); }

	//	if (!atoms.entries.empty()) {
	//		file << atoms.composeString();
	//		file << singlebonds.composeString();
	//		file << pairs.composeString();
	//		file << anglebonds.composeString();
	//		file << dihedralbonds.composeString();
	//		file << improperdihedralbonds.composeString();
	//	}

	//	if (path.extension() == ".top")
	//		file << "[ system ]\n" << system << "\n";

	//	//if (!molecules.entries.empty()) { file << molecules.composeString(); }
	//	// If we have the same submolecule multiple times in a row, we only print it once together with a count of how many there are
	//	if (!molecules.entries.empty())
	//		file << molecules.title << "\n" << molecules.legend << "\n";
	//	for (int i = 0; i < molecules.entries.size(); i++) {
	//		std::ostringstream oss;
	//		int count = 1;
	//		while (i + 1 < molecules.entries.size() && molecules.entries[i].name == molecules.entries[i + 1].name) {
	//			count++;
	//			i++;
	//		}
	//		file << molecules.entries[i].includeTopologyFile->name << " " << count << "\n";
	//	}
	//}

	//// Also cache the file
	//WriteFileToBinaryCache(*this, path);

	//for (const auto& [name, include] : includeTopologies) {
	//	include->printToFile(path.parent_path() / ("topol_" + name + ".itp"), printForcefieldinclude);
	//}
}
void TopologyFile::AppendTopology(const std::shared_ptr<TopologyFile>& other) {

	//if (!other->moleculetype.has_value()) {
	//	throw std::runtime_error(std::format("Trying to append a topology file without a moleculetype section: {}", other->path.string()));
	//}

	//if (!includeTopologies.contains(other->moleculetype.value().name))
	//	includeTopologies.insert({ other->moleculetype.value().name, other });

	//const int globalIndexOfFirstParticle = molecules.entries.empty()
	//	? 0 : molecules.entries.back().GlobalIndexOfFinalParticle() + 1;
	//molecules.entries.emplace_back(other->name, includeTopologies.at(other->name), globalIndexOfFirstParticle);
}
//std::optional<TopologyFile::ForcefieldInclude> TopologyFile::ThisAndAllSubmoleculesShareTheSameForcefieldinclude() const {
//
//	//std::optional<TopologyFile::ForcefieldInclude> sharedff = std::nullopt;
//
//	//for (const auto& mol : GetAllSubMolecules()) {
//	//	if (!mol.includeTopologyFile->forcefieldInclude)
//	//		continue;
//
//	//	if (!sharedff)
//	//		sharedff = mol.includeTopologyFile->forcefieldInclude;
//	//	else if (sharedff.value().name != mol.includeTopologyFile->forcefieldInclude.value().name)
//	//		return std::nullopt;
//	//}
//	//return sharedff;
//}

//
//TopologyFile::SectionRange<TopologyFile::AtomsEntry> TopologyFile::GetAllElements<TopologyFile::AtomsEntry>() const {
//	return TopologyFile::SectionRange<AtomsEntry>(this);
//}
//TopologyFile::SectionRange<TopologyFile::SingleBond> TopologyFile::GetAllSinglebonds() const {
//	return TopologyFile::SectionRange<SingleBond>(this);
//}
//TopologyFile::SectionRange<TopologyFile::Pair> TopologyFile::GetAllPairs() const {
//	return TopologyFile::SectionRange<Pair>(this);
//}
//TopologyFile::SectionRange<TopologyFile::AngleBond> TopologyFile::GetAllAnglebonds() const {
//	return TopologyFile::SectionRange<AngleBond>(this);
//}
//TopologyFile::SectionRange<TopologyFile::DihedralBond> TopologyFile::GetAllDihedralbonds() const {
//	return TopologyFile::SectionRange<DihedralBond>(this);
//}
//TopologyFile::SectionRange<TopologyFile::ImproperDihedralBond> TopologyFile::GetAllImproperDihedralbonds() const {
//	return TopologyFile::SectionRange<ImproperDihedralBond>(this);
//}
//TopologyFile::SectionRange<TopologyFile::MoleculeEntry> TopologyFile::GetAllSubMolecules() const {
//	return TopologyFile::SectionRange<MoleculeEntry>(this);
//}


void MDFiles::MergeFiles(GroFile& leftGro, TopologyFile& leftTop, GroFile& rightGro, std::shared_ptr<TopologyFile> rightTop) 
{
	int newGroId = leftGro.atoms.back().gro_id + 1;
	// First merge the GRO, and get a mapping to the new IDs
	for (const GroRecord& atom : rightGro.atoms) {
		leftGro.atoms.emplace_back(atom);
		leftGro.atoms.back().gro_id = newGroId % 100'000;
		newGroId++;
	}

	// To merge the top files, we simply need to add it as an include topology
	leftTop.AppendTopology(rightTop);	
}

SimulationFilesCollection::SimulationFilesCollection(const fs::path& workDir) {
	grofile = std::make_unique<GroFile>(workDir / "molecule/conf.gro");
	topfile = std::make_unique<TopologyFile>(workDir / "molecule/topol.top");
}


PDBfile::PDBfile(const fs::path& path) : mPath(path) {
	if (!(path.extension().string() == std::string{ ".pdb" }))
		throw std::runtime_error("Expected .pdb extension");
	if (!fs::exists(path))
		throw std::runtime_error(std::format("File \"{}\" was not found", path.string()));



	std::ifstream file;
	file.open(path);
	if (!file.is_open() || file.fail()) {
		throw std::runtime_error(std::format("Failed to open file {}\n", path.string()).c_str());
	}

	std::string line{};
	while (getline(file, line)) {
		if (line.length() < 80) continue;

		bool isATOM = (line.substr(0, 4) == "ATOM");
		bool isHETATM = (line.substr(0, 6) == "HETATM");

		if (isATOM || isHETATM) {
			ATOM atom;

			atom.atomSerialNumber = std::stoi(line.substr(6, 5));
			strncpy(atom.atomName, line.substr(12, 4).c_str(), 4);
			atom.altLocIndicator = line[16];
			strncpy(atom.resName, line.substr(17, 3).c_str(), 3);
			atom.chainID = line[21];
			atom.resSeq = std::stoi(line.substr(22, 4));
			atom.iCode = line[26];

			// Convert coordinates from Ångströms to nanometers
			atom.position.x = std::stof(line.substr(30, 8)) * 0.1f;
			atom.position.y = std::stof(line.substr(38, 8)) * 0.1f;
			atom.position.z = std::stof(line.substr(46, 8)) * 0.1f;

			atom.occupancy = std::stof(line.substr(54, 6));
			atom.tempFactor = std::stof(line.substr(60, 6));
			strncpy(atom.segmentIdentifier, line.substr(72, 4).c_str(), 4);
			atom.elementSymbol = line[76];
			strncpy(atom.charge, line.substr(78, 2).c_str(), 2);

			if (isATOM) {
				ATOMS.push_back(atom);
			}
			else if (isHETATM) {
				HETATMS.push_back(atom);
			}
		}
	}
	
}
