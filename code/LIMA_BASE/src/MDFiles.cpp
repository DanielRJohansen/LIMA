#include "MDFiles.h"
#include "Filehandling.h"
#include "MDFilesSerialization.h"

#include <algorithm>
#include <format>


using namespace Filehandler;
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

void TopologyFile::MoleculeEntry::composeString(std::ostringstream& oss) const {
	if (includeTopologyFile->name == "")
		throw std::runtime_error("Trying to save a topology with an empty include topology name");
	oss << includeTopologyFile->name;
}

void TopologyFile::MoleculetypeEntry::composeString(std::ostringstream& oss) const {
	oss << std::right << std::setw(10) << name << std::setw(10) << nrexcl;
}

void TopologyFile::AtomsEntry::composeString(std::ostringstream& oss) const {
	if (section_name) {
		oss << section_name.value() << "\n";
	}
	oss << std::right
		<< std::setw(10) << id
		<< std::setw(10) << type
		<< std::setw(10) << resnr
		<< std::setw(10) << residue
		<< std::setw(10) << atomname
		<< std::setw(10) << cgnr
		<< std::setw(10) << std::fixed << std::setprecision(2) << charge
		<< std::setw(10) << std::fixed << std::setprecision(3) << mass;
}


void LoadDefaultIncludeTopologies(std::unordered_map<std::string, LazyLoadFile<TopologyFile>>& includedFiles, fs::path srcDir) {
	if (!fs::exists(srcDir) || !fs::is_directory(srcDir)) {
        throw std::invalid_argument("Source path is not a directory: " + srcDir.string());
    }
	for (const auto& entry : fs::directory_iterator(srcDir)) {
		if (entry.is_directory()) {
			fs::path subDir = entry.path();
			fs::path itpFile = subDir / (subDir.filename().string() + ".itp");
			if (fs::exists(itpFile)) {
				includedFiles.insert({ subDir.filename().string(), {itpFile} });
			}
		}
	}
}



bool isFloat(const std::string& str) {
	std::istringstream iss(str);
	float f;
	iss >> std::noskipws >> f; // noskipws considers leading/trailing whitespace invalid
	return iss.eof() && !iss.fail();
}

GroFile::GroFile(const fs::path& path) : m_path(path){
	if (!(path.extension().string() == std::string{ ".gro" }))
		throw std::runtime_error("Expected .gro extension");
	if (!fs::exists(path))
		throw std::runtime_error(std::format("File \"{}\" was not found", path.string()));

	lastModificationTimestamp = fs::last_write_time(path);

	if (UseCachedBinaryFile(path, lastModificationTimestamp)) {
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
			}
			prevLine = line;
		}

		if (!prevLine.empty()) {
			int dim = 0;
			std::stringstream ss(prevLine);
			while (std::getline(ss, word, ' ')) {
				if (!word.empty()) {
					box_size[dim++] = std::stof(word);
				}
			}
		}

		assert(atoms.size() == nAtoms);

		// Save a binary cached version of the file to we can read it faster next time
		WriteFileToBinaryCache(*this);
	}
}


class TopologySectionGetter {
	int dihedralCount = 0;
	int dihedraltypesCount = 0;
public:
	TopologySection operator()(const std::string& directive) {
		if (directive == "molecules") return molecules;
		if (directive == "moleculetype") return moleculetype;
		if (directive == "atoms") return atoms;
		if (directive == "bonds") return bonds;
		if (directive == "pairs") return pairs;
		if (directive == "angles") return angles;
		if (directive == "dihedrals") {
			if (++dihedralCount <= 2) return dihedralCount == 1 ? dihedrals : impropers;
			throw std::runtime_error("Encountered 'dihedral' directive more than 2 times in .itp/.top file");
		}
		if (directive == "position_restraints") return position_restraints;
		if (directive == "system") return _system;
		if (directive == "cmap") return cmap;

		if (directive == "atomtypes") return atomtypes;
		if (directive == "pairtypes") return pairtypes;
		if (directive == "bondtypes") return bondtypes;
		if (directive == "constrainttypes") return constainttypes;
		if (directive == "angletypes") return angletypes;
		if (directive == "dihedraltypes") {
			if (++dihedraltypesCount <= 2) return dihedraltypesCount == 1 ? dihedraltypes : impropertypes;
			throw std::runtime_error("Encountered 'dihedraltypes' directive more than 2 times in .itp/ file");
		}
		if (directive == "impropertypes") return impropertypes;

		throw std::runtime_error(std::format("Got unexpected topology directive: {}", directive));
	}
};

std::string extractSectionName(const std::string& line) {
	size_t start = line.find('[');
	size_t end = line.find(']', start);
	if (start != std::string::npos && end != std::string::npos) {
		// Extract the text between '[' and ']'
		std::string sectionName = line.substr(start + 1, end - start - 1);
		removeWhitespace(sectionName);
		return sectionName;
	}
	return "";
}


std::string GetCleanFilename(const fs::path& path) {
	auto filename = path.stem().string();
	const std::string prefix = "topol_";
	return filename.starts_with(prefix) ? filename.substr(prefix.length()) : filename;
}

fs::path SearchForFile(fs::path dir, const std::string& filename) {
	const std::array<std::string, 2> extensions = { ".itp", ".top" };
	const std::array<std::string, 2> prefixes = { std::string(""), "topol_" };

	for (const auto& ext : extensions) {
		for (const auto& prefix : prefixes) {
			fs::path path = dir / (prefix + filename + ext);
			if (fs::exists(path)) {
				return path;
			}
		}
	}

	throw std::runtime_error(std::format("Could not find file \"{}\" in directory \"{}\"", filename, dir.string()));
}

TopologyFile::TopologyFile() {
	LoadDefaultIncludeTopologies(includedFiles, Filehandler::GetLimaDir() / "resources/Lipids");
}

template <int n>
bool VerifyAllParticlesInBondExists(const std::vector<int>& groIdToLimaId, int ids[n]) {
	for (int i = 0; i < n; i++) {
		if (ids[i] >= groIdToLimaId.size() || groIdToLimaId[ids[i]] == -1)
			return false;
	}
	return true;
}

/// <summary></summary>
/// <param name="line"></param>
/// <param name="currentSection"></param>
/// <param name="sectionGetter"></param>
/// <returns>True if we change section/stay on no section, in which case we should 'continue' to the next line. False otherwise</returns>
bool HandleTopologySectionStartAndStop(const std::string& line, TopologySection& currentSection, TopologySectionGetter& sectionGetter) {
	if (line.size() == 0) {
		currentSection = no_section;
		return true;
	}
	else if (std::all_of(line.begin(), line.end(), isspace)) {	// Check if all chars in line is space
		return true;	// TODO Set to nosection here
	}
	else if (line[0] == '[') {
		currentSection = sectionGetter(extractSectionName(line));
		return true;
	}
	
	return false; // We did not change section, and should continue parsing the current line
}

TopologyFile::TopologyFile(const fs::path& path) : path(path), name(GetCleanFilename(path))
{
	LoadDefaultIncludeTopologies(includedFiles, Filehandler::GetLimaDir() / "resources/Lipids");
	if (!(path.extension().string() == std::string{ ".top" } || path.extension().string() == ".itp"))
		throw std::runtime_error("Expected .top or .itp extension");
	if (!fs::exists(path))
		throw std::runtime_error(std::format("File \"{}\" was not found", path.string()));

	lastModificationTimestamp = fs::last_write_time(path);

	if (UseCachedBinaryFile(path, lastModificationTimestamp) && ENABLE_FILE_CACHING) {
		readTopFileFromBinaryCache(path, *this);

		// Any include topologies will not be in this particular cached binary, so iterate through them, and load them
		for (auto& molecule : molecules.entries) {
			assert(molecule.includeTopologyFile == nullptr);
			assert(molecule.name != "");

			if (includedFiles.count(molecule.name) == 0)
				includedFiles.emplace(molecule.name, SearchForFile(path.parent_path(), molecule.name));
			molecule.includeTopologyFile = includedFiles.at(molecule.name).Get();
		}
	}
	else {
		std::ifstream file;
		file.open(path);
		if (!file.is_open() || file.fail()) {
			throw std::runtime_error(std::format("Failed to open file {}\n", path.string()));
		}

		TopologySection current_section{ TopologySection::title };
		TopologySectionGetter getTopolSection{};

		std::string line{}, word{};
		std::string sectionname = "";
		std::vector<int> groIdToLimaId;


		while (getline(file, line)) {
			if (HandleTopologySectionStartAndStop(line, current_section, getTopolSection)) {
				continue;
			}

			// Check if current line is commented
			if (firstNonspaceCharIs(line, commentChar) && current_section != TopologySection::title && current_section != TopologySection::atoms) { continue; }	// Only title-sections + atoms reads the comments

			std::istringstream iss(line);
			switch (current_section)
			{
			case TopologySection::title:
				title.append(line + "\n");	// +\n because getline implicitly strips it away.
				break;
			case TopologySection::molecules: {
				// TODO: This is incorrect, we should have a section where we look for "#include" and include those files
				std::string include_name;
				iss >> include_name;

				// Handle the case where there simply a random SOL that does not refer to a file.. Terrible standard...
				if (include_name == "SOL")
					continue;


				if (includedFiles.count(include_name) == 0) {
					const fs::path includePath = SearchForFile(path.parent_path(), include_name);
					includedFiles.emplace(include_name, includePath);
				}

				molecules.entries.emplace_back(include_name, includedFiles.at(include_name).Get());
			}
			case moleculetype:
			{
				TopologyFile::MoleculetypeEntry moleculetype{};
				iss >> moleculetype.name >> moleculetype.nrexcl;
				moleculetypes.entries.emplace_back(moleculetype);
				break;
			}
			case TopologySection::atoms:
			{
				if (firstNonspaceCharIs(line, ';')) {
					// TODO: Test for residue or lipid_section in the [1] position of the comment instead
					if (iss.str().find("type") == std::string::npos) {
						sectionname = line;
					}
				}
				else {
					TopologyFile::AtomsEntry atom;
					int groId;
					iss >> groId>> atom.type >> atom.resnr >> atom.residue >> atom.atomname >> atom.cgnr >> atom.charge >> atom.mass;

					if (groIdToLimaId.size() < groId + 1)
						groIdToLimaId.resize(groId + 1, -1);
					groIdToLimaId[groId] = atoms.entries.size();
					atom.id = groIdToLimaId[groId];
					atoms.entries.emplace_back(atom);

					if (sectionname != "") {
						atoms.entries.back().section_name = sectionname;
						sectionname = "";
					}
				}
				break;
			}
			case TopologySection::bonds: {
				TopologyFile::SingleBond singlebond{};
				int groIds[2];
				iss >> groIds[0] >> groIds[1] >> singlebond.funct;
				if (!VerifyAllParticlesInBondExists<2>(groIdToLimaId, groIds))
					break;
				for (int i = 0; i < 2; i++) 
					singlebond.ids[i] = groIdToLimaId[groIds[i]];
				if (singlebond.ids[0] == 1770 && singlebond.ids[1] == 1772) {
					auto a = name;
					int c = 0;
				}
				singlebonds.entries.emplace_back(singlebond);
				break;
			}
			case TopologySection::pairs: {
				TopologyFile::Pair pair{};
				int groIds[2];
				iss >> groIds[0] >> groIds[1] >> pair.funct;
				if (!VerifyAllParticlesInBondExists<2>(groIdToLimaId, groIds))
					break;
				for (int i = 0; i < 2; i++)
					pair.ids[i] = groIdToLimaId.at(groIds[i]);
				pairs.entries.emplace_back(pair);
				break;
			}
			case TopologySection::angles: {
				TopologyFile::AngleBond angle{};
				int groIds[3];
				iss >> groIds[0] >> groIds[1] >> groIds[2] >> angle.funct;
				if (!VerifyAllParticlesInBondExists<3>(groIdToLimaId, groIds))
					break;
				for (int i = 0; i < 3; i++)					
					angle.ids[i] = groIdToLimaId.at(groIds[i]);								
				anglebonds.entries.emplace_back(angle);
				break;
			}
			case TopologySection::dihedrals: {
				TopologyFile::DihedralBond dihedral{};
				int groIds[4];
				iss >> groIds[0] >> groIds[1] >> groIds[2] >> groIds[3] >> dihedral.funct;
				if (!VerifyAllParticlesInBondExists<4>(groIdToLimaId, groIds))
					break;
				for (int i = 0; i < 4; i++)					
					dihedral.ids[i] = groIdToLimaId.at(groIds[i]);				
				dihedralbonds.entries.emplace_back(dihedral);
				break;
			}
			case TopologySection::impropers: {
				TopologyFile::ImproperDihedralBond improper{};
				int groIds[4];
				iss >> groIds[0] >> groIds[1] >> groIds[2] >> groIds[3] >> improper.funct;
				if (!VerifyAllParticlesInBondExists<4>(groIdToLimaId, groIds))
					break;
				for (int i = 0; i < 4; i++)
					improper.ids[i] = groIdToLimaId.at(groIds[i]);								
				improperdihedralbonds.entries.emplace_back(improper);
				improperdihedralbonds.entries.back().sourceLine = line;
				break;
			}
			default:
				// Do nothing
				//throw std::runtime_error("Illegal state");
				break;
			}
		}
		WriteFileToBinaryCache(*this);
	}


	//// Verify that atoms id's are a sequence starting at 1
	//for (int i = 0; i < atoms.entries.size(); i++) {
	//	if (atoms.entries[i].nr != i+1) {
	//		throw std::runtime_error("Atoms are not in sequence starting at 1");
	//	}
	//}
	//for (const auto& mol : molecules.entries) {
	//	for (int i = 0; i < mol.includeTopologyFile->atoms.entries.size(); i++) {
	//		if (mol.includeTopologyFile->atoms.entries[i].nr != i+1) {
	//			throw std::runtime_error("Atoms are not in sequence starting at 1");
	//		}
	//	}
	//}
}

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








std::string TopologyFile::generateLegend(std::vector<std::string> elements)
{
	std::ostringstream legend;
	legend << ';'; // Start with a semicolon

	for (const auto& element : elements) {
		legend << std::setw(10) << std::right << element;
	}
	return legend.str();
}

void TopologyFile::printToFile(const std::filesystem::path& path) const {
	const auto ext = path.extension().string();
	if (ext != ".top" && ext != ".itp") { throw std::runtime_error(std::format("Got {} extension, expectec [.top/.itp]", ext)); }

	std::ofstream file(path);
	if (!file.is_open()) {
		throw std::runtime_error(std::format("Failed to open file {}", path.string()));
	}

	file << title << "\n\n";

	// TODO: What i currently use "molecules" for should actually be an "#include" argument..

	if (!molecules.entries.empty()) { file << molecules.composeString(); }
	if (!moleculetypes.entries.empty()) { file << moleculetypes.composeString(); }

	
	file << atoms.composeString();
	file << singlebonds.composeString();
	file << pairs.composeString();
	file << anglebonds.composeString();
	file << dihedralbonds.composeString();
	file << improperdihedralbonds.composeString();



	// Also cache the file
	WriteFileToBinaryCache(*this, path);
}


TopologyFile::SectionRange<TopologyFile::AtomsEntry> TopologyFile::GetAllAtoms() const {
	return TopologyFile::SectionRange<AtomsEntry>(this);
}
TopologyFile::SectionRange<TopologyFile::SingleBond> TopologyFile::GetAllSinglebonds() const {
	return TopologyFile::SectionRange<SingleBond>(this);
}
TopologyFile::SectionRange<TopologyFile::Pair> TopologyFile::GetAllPairs() const {
	return TopologyFile::SectionRange<Pair>(this);
}
TopologyFile::SectionRange<TopologyFile::AngleBond> TopologyFile::GetAllAnglebonds() const {
	return TopologyFile::SectionRange<AngleBond>(this);
}
TopologyFile::SectionRange<TopologyFile::DihedralBond> TopologyFile::GetAllDihedralbonds() const {
	return TopologyFile::SectionRange<DihedralBond>(this);
}
TopologyFile::SectionRange<TopologyFile::ImproperDihedralBond> TopologyFile::GetAllImproperDihedralbonds() const {
	return TopologyFile::SectionRange<ImproperDihedralBond>(this);
}
TopologyFile::SectionRange<TopologyFile::MoleculeEntry> TopologyFile::GetAllSubMolecules() const {
	return TopologyFile::SectionRange<MoleculeEntry>(this);
}


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

GenericItpFile::GenericItpFile(const fs::path& path) {
	if (path.extension().string() != ".itp") { throw std::runtime_error(std::format("Expected .itp extension with file {}", path.string())); }
	if (!fs::exists(path)) { throw std::runtime_error(std::format("File \"{}\" was not found", path.string())); }

	std::ifstream file;
	file.open(path);
	if (!file.is_open() || file.fail()) {
		throw std::runtime_error(std::format("Failed to open file {}\n", path.string()));
	}

	//lastModificationTimestamp = fs::last_write_time(path);

	//if (UseCachedBinaryFile(path, lastModificationTimestamp) && ENABLE_FILE_CACHING) {
	//	readTopFileFromBinaryCache(path, *this);

	//	// Any include topologies will not be in this particular cached binary, so iterate through them, and load them
	//	for (auto& molecule : molecules.entries) {
	//		assert(molecule.includeTopologyFile == nullptr);
	//		assert(molecule.name != "");

	//		if (includedFiles.count(molecule.name) == 0)
	//			includedFiles.emplace(molecule.name, SearchForFile(path.parent_path(), molecule.name));
	//		molecule.includeTopologyFile = includedFiles.at(molecule.name).Get();
	//	}
	//}
	//else {


	TopologySection current_section{ TopologySection::title };
	TopologySectionGetter getTopolSection{};
	bool newSection = true;

	std::string line{}, word{};
	std::string sectionname = "";
	std::vector<int> groIdToLimaId;


	while (getline(file, line)) {
		if (HandleTopologySectionStartAndStop(line, current_section, getTopolSection)) {
			newSection = true;
			continue;
		}

		// Check if current line is commented
		if (firstNonspaceCharIs(line, TopologyFile::commentChar) && current_section != TopologySection::title && current_section != TopologySection::atoms) { continue; }	// Only title-sections + atoms reads the comments

		// Check if this line contains another illegal keyword
		if (firstNonspaceCharIs(line, '#')) {// Skip ifdef, include etc TODO: this needs to be implemented at some point
			if (line.find("#include") != std::string::npos) {
				throw std::runtime_error("This class should never be reading files containing thye #include keyword!");
			}
			continue;
		}


		if (newSection) {
			if (GetSection(current_section).has_value()) {
				throw std::runtime_error("Found the same section muliple times in the same file");
			}
			sections.insert({ current_section, {} });
			newSection = false;
		}

		GetSection(current_section)->get().lines.emplace_back(line);	// OPTIM: i prolly should cache this address instead
	}
}