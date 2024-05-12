#include "MDFiles.h"
#include "Filehandling.h"
#include "MDFilesSerialization.h"

#include <algorithm>
#include <format>
#include <fstream>
#include <future>

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
	record.residue_name = trimSpaces(line.substr(5, 5));

	// Parse atom name (5 characters)
	record.atom_name = trimSpaces(line.substr(10, 5));

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
		<< std::setw(5) << std::left << record.residue_name
		<< std::setw(5) << std::left << record.atom_name
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

void ParsedTopologyFile::MoleculeEntry::composeString(std::ostringstream& oss) const {
	if (includeTopologyFile->name == "")
		throw std::runtime_error("Trying to save a topology with an empty include topology name");
	oss << includeTopologyFile->name;
}

void ParsedTopologyFile::MoleculetypeEntry::composeString(std::ostringstream& oss) const {
	oss << std::right << std::setw(10) << name << std::setw(10) << nrexcl;
}

void ParsedTopologyFile::AtomsEntry::composeString(std::ostringstream& oss) const {
	if (section_name) {
		oss << section_name.value() << "\n";
	}
	oss << std::right
		<< std::setw(10) << nr
		<< std::setw(10) << type
		<< std::setw(10) << resnr
		<< std::setw(10) << residue
		<< std::setw(10) << atomname
		<< std::setw(10) << cgnr
		<< std::setw(10) << std::fixed << std::setprecision(2) << charge
		<< std::setw(10) << std::fixed << std::setprecision(3) << mass;
}


void LoadDefaultIncludeTopologies(std::unordered_map<std::string, LazyLoadFile<ParsedTopologyFile>>& includedFiles, fs::path srcDir) {
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

ParsedGroFile::ParsedGroFile(const fs::path& path) : m_path(path){
	if (!(path.extension().string() == std::string{ ".gro" }))
		throw std::runtime_error("Expected .gro extension");
	if (!fs::exists(path))
		throw std::runtime_error(std::format("File \"{}\" was not found", path.string()));

	lastModificationTimestamp = fs::last_write_time(path);

	if (UseCachedBinaryFile(path, lastModificationTimestamp)) {
		*this = readGroFileFromBinaryCache(path);
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
					*box_size.placeAt(dim++) = std::stof(word);
				}
			}
		}

		assert(atoms.size() == nAtoms);

		// Save a binary cached version of the file to we can read it faster next time
		WriteFileToBinaryCache(*this);
	}
}

enum TopologySection {title, molecules, moleculetype, atoms, bonds, pairs, angles, dihedrals, impropers, position_restraints, _system, cmap, no_section};

class TopologySectionGetter {
	int dihedralCount = 0;

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

ParsedTopologyFile::ParsedTopologyFile() {
	LoadDefaultIncludeTopologies(includedFiles, Filehandler::GetLimaDir() / "resources/Lipids");
}

ParsedTopologyFile::ParsedTopologyFile(const fs::path& path) : path(path), name(GetCleanFilename(path))
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
			//includeTopology.includeTopologyFile = std::make_shared<ParsedTopologyFile>(includeTopology.name);
			// 
			// 
			if (includedFiles.count(molecule.name) == 0)
				includedFiles.emplace(molecule.name, SearchForFile(path.parent_path(), molecule.name));
			molecule.includeTopologyFile = includedFiles.at(molecule.name).Get();
			//includeTopology.includeTopologyFile = includedFiles.at(includeTopology.name).Get().get();
		}
		//LoadAllSubMolecules();

		return;
	}



	assert(path.extension().string() == std::string{ ".top" } || path.extension().string() == ".itp");

	const char delimiter = ' ';
	const std::vector<char> ignores = { ';', '#' };

	std::ifstream file;
	file.open(path);
	if (!file.is_open() || file.fail()) {
		throw std::runtime_error(std::format("Failed to open file {}\n", path.string()));
	}

	TopologySection current_section{ TopologySection::title };
	TopologySectionGetter getTopolSection{};

	std::string line{}, word{};
	std::string sectionname = "";
//	std::vector<std::future<std::pair<std::string, std::unique_ptr<ParsedTopologyFile>>>> includeTopologies;	// <name, file>

	while (getline(file, line)) {

		if (line.size() == 0) {
			current_section = no_section;
			continue;
		}

		if (line[0] == '[') {
			current_section = getTopolSection(extractSectionName(line));
			continue;
		}

		// Check if current line is commented
		if (firstNonspaceCharIs(line, ';') && current_section != TopologySection::title && current_section != TopologySection::atoms) { continue; }	// Only title-sections + atoms reads the comments

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
				//const fs::path include_path = fs::exists(path.parent_path() / (include_name + ".itp"))
				//	? path.parent_path() / (include_name + ".itp")
				//	: path.parent_path() / ("topol_" + include_name + ".itp");


				const fs::path includePath = SearchForFile(path.parent_path(), include_name);
				includedFiles.emplace(include_name, includePath);
			}

			molecules.entries.emplace_back(include_name, includedFiles.at(include_name).Get());
		}
		case moleculetype:
		{
			ParsedTopologyFile::MoleculetypeEntry moleculetype{};
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
				ParsedTopologyFile::AtomsEntry atom;
				iss >> atom.nr >> atom.type >> atom.resnr >> atom.residue >> atom.atomname >> atom.cgnr >> atom.charge >> atom.mass;
				atoms.entries.emplace_back(atom);

				if (sectionname != "") {
					atoms.entries.back().section_name = sectionname;
					sectionname = "";
				}
			}
			break;
		}
		case TopologySection::bonds: {
			ParsedTopologyFile::SingleBond singlebond{};
			iss >> singlebond.atomGroIds[0] >> singlebond.atomGroIds[1] >> singlebond.funct;
			singlebonds.entries.emplace_back(singlebond);
			break;
		}
		case TopologySection::pairs: {
			ParsedTopologyFile::Pair pair{};
			iss >> pair.atomGroIds[0] >> pair.atomGroIds[1] >> pair.funct;
			pairs.entries.emplace_back(pair);
			break;
		}
		case TopologySection::angles: {
			ParsedTopologyFile::AngleBond angle{};
			iss >> angle.atomGroIds[0] >> angle.atomGroIds[1] >> angle.atomGroIds[2] >> angle.funct;
			anglebonds.entries.emplace_back(angle);
			break;
		}
		case TopologySection::dihedrals: {
			ParsedTopologyFile::DihedralBond dihedral{};
			iss >> dihedral.atomGroIds[0] >> dihedral.atomGroIds[1] >> dihedral.atomGroIds[2] >> dihedral.atomGroIds[3] >> dihedral.funct;
			dihedralbonds.entries.emplace_back(dihedral);
			break;
		}
		case TopologySection::impropers: {
			ParsedTopologyFile::ImproperDihedralBond improper{};
			iss >> improper.atomGroIds[0] >> improper.atomGroIds[1] >> improper.atomGroIds[2] >> improper.atomGroIds[3] >> improper.funct;
			improperdihedralbonds.entries.emplace_back(improper);
			break;
		}
		default:
			// Do nothing
			//throw std::runtime_error("Illegal state");
			break;
		}
	}

	//LoadAllSubMolecules();

	WriteFileToBinaryCache(*this);
}

void ParsedGroFile::printToFile(const std::filesystem::path& path) const {
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

//void ParsedGroFile::addEntry(std::string residue_name, std::string atom_name, const Float3& position) {
//	if (atoms.size() == 0) {
//		atoms.push_back({ 0, residue_name, atom_name, 0, position, std::nullopt });
//	}
//	else {
//		atoms.push_back({atoms.back().residue_number+1, residue_name, atom_name, atoms.back().gro_id+1, position, std::nullopt});
//	}
//}










std::string ParsedTopologyFile::generateLegend(std::vector<std::string> elements)
{
	std::ostringstream legend;
	legend << ';'; // Start with a semicolon

	for (const auto& element : elements) {
		legend << std::setw(10) << std::right << element;
	}
	return legend.str();
}

void ParsedTopologyFile::printToFile(const std::filesystem::path& path) const {
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


//void ParsedTopologyFile::IncrementIds(int atomNrIncrement, int resNrIncrement) {
//
//	for (auto& atom : atoms.entries) {
//		atom.nr += atomNrIncrement;
//		atom.resnr += resNrIncrement;
//	}
//	for (auto& singlebond : singlebonds.entries) {
//		singlebond.IncrementIds(atomNrIncrement);
//	}
//	for (auto& pair : pairs.entries) {
//		pair.IncrementIds(atomNrIncrement);
//	}
//	for (auto& anglebond : anglebonds.entries) {
//		anglebond.IncrementIds(atomNrIncrement);
//	}
//	for (auto& dihedralbond : dihedralbonds.entries) {
//		dihedralbond.IncrementIds(atomNrIncrement);
//	}
//	for (auto& improperdihedralbond : improperdihedralbonds.entries) {
//		improperdihedralbond.IncrementIds(atomNrIncrement);
//	}
//
//	for (auto& molecule : molecules.entries) {
//		molecule.includeTopologyFile->IncrementIds(atomNrIncrement, resNrIncrement);
//	}
//}

ParsedTopologyFile::SectionRange<ParsedTopologyFile::AtomsEntry> ParsedTopologyFile::GetAllAtoms() const {
	return ParsedTopologyFile::SectionRange<AtomsEntry>(this);
}
//ParsedTopologyFile::SectionRange<ParsedTopologyFile::SingleBond> ParsedTopologyFile::GetSinglebonds() {
//	return ParsedTopologyFile::SectionRange<SingleBond>(this);
//}
ParsedTopologyFile::SectionRange<ParsedTopologyFile::SingleBond> ParsedTopologyFile::GetAllSinglebonds() const {
	return ParsedTopologyFile::SectionRange<SingleBond>(this);
}
ParsedTopologyFile::SectionRange<ParsedTopologyFile::Pair> ParsedTopologyFile::GetAllPairs() const {
	return ParsedTopologyFile::SectionRange<Pair>(this);
}
ParsedTopologyFile::SectionRange<ParsedTopologyFile::AngleBond> ParsedTopologyFile::GetAllAnglebonds() const {
	return ParsedTopologyFile::SectionRange<AngleBond>(this);
}
ParsedTopologyFile::SectionRange<ParsedTopologyFile::DihedralBond> ParsedTopologyFile::GetAllDihedralbonds() const {
	return ParsedTopologyFile::SectionRange<DihedralBond>(this);
}
ParsedTopologyFile::SectionRange<ParsedTopologyFile::ImproperDihedralBond> ParsedTopologyFile::GetAllImproperDihedralbonds() const {
	return ParsedTopologyFile::SectionRange<ImproperDihedralBond>(this);
}
ParsedTopologyFile::SectionRange<ParsedTopologyFile::MoleculeEntry> ParsedTopologyFile::GetAllSubMolecules() const {
	return ParsedTopologyFile::SectionRange<MoleculeEntry>(this);
}


void MDFiles::MergeFiles(ParsedGroFile& leftGro, ParsedTopologyFile& leftTop, ParsedGroFile& rightGro, std::unique_ptr<ParsedTopologyFile> rightTop) 
{
	// First merge the GRO, and get a mapping to the new IDs
	for (const GroRecord& atom : rightGro.atoms) {
		leftGro.atoms.emplace_back(atom);
	}



	// To merge the top files, we simply need to add it as an include topology
	leftTop.AppendTopology(std::move(rightTop));	
}


















ParsedLffFile::LffSection getLffSection(const std::string& directive) {
	if (directive == "singlebonds")
		return ParsedLffFile::singlebond;
	else if (directive == "anglebonds")
		return ParsedLffFile::anglebond;
	else if (directive == "dihedralbonds")
		return ParsedLffFile::dihedralbond;
	else if (directive == "improperdihedralbonds")
		return ParsedLffFile::improperdihedralbond;
	else
		return ParsedLffFile::no_section;
}

ParsedLffFile::ParsedLffFile(const fs::path& path) : path(path) {
	if (path.extension().string() != std::string{ ".lff" }) { throw std::runtime_error(std::format("Expected .lff extension with file {}", path.string())); }
	if (!fs::exists(path)) { throw std::runtime_error(std::format("File \"{}\" was not found", path.string())); }


	std::ifstream file;
	file.open(path);
	if (!file.is_open() || file.fail()) {
		throw std::runtime_error(std::format("Failed to open file {}\n", path.string()));
	}

	LffSection current_section{ title };

	// Forward declaring for optimization reasons
	std::string line{}, word{};
	std::string sectionname = "";
	while (getline(file, line)) {

		if (line.size() == 0) {
			current_section = no_section;
			continue;
		}

		if (line[0] == '#') {
			//current_section = getNextTopologySection(current_section, extractSectionName(line));
			current_section = getLffSection(line.substr(2));
			continue;
		}

		// Check if current line is commented
		if (firstNonspaceCharIs(line, '/') && current_section != title) { continue; }	// Only title-sections reads the comments


		std::istringstream iss(line);
		switch (current_section)
		{
		case singlebond: {
			ParsedLffFile::Singlebond singlebond{};
			std::string atomtypes[2];
			iss >> singlebond.global_ids[0] >> singlebond.global_ids[1] >> atomtypes[0] >> atomtypes[1] >> singlebond.b0 >> singlebond.kb;
			singlebonds.entries.emplace_back(singlebond);
			break;
		}
		case anglebond: {
			ParsedLffFile::Anglebond anglebond{};
			std::string atomtypes[3];
			iss >> anglebond.global_ids[0] >> anglebond.global_ids[1] >> anglebond.global_ids[2] >> atomtypes[0] >> atomtypes[1] >> atomtypes[2] >> anglebond.theta0 >> anglebond.ktheta;
			anglebonds.entries.emplace_back(anglebond);
			break;
		}
		case dihedralbond: {
			ParsedLffFile::Dihedralbond dihedralbond{};
			std::string atomtypes[4];
			iss >> dihedralbond.global_ids[0] >> dihedralbond.global_ids[1] >> dihedralbond.global_ids[2] >> dihedralbond.global_ids[3] >> atomtypes[0] >> atomtypes[1] >> atomtypes[2] >> atomtypes[3] >> dihedralbond.phi0 >> dihedralbond.kphi >> dihedralbond.n;
			dihedralbonds.entries.emplace_back(dihedralbond);
			break;
		}
		case improperdihedralbond: {
			ParsedLffFile::Improperdihedralbond improperdihedralbond{};
			std::string atomtypes[4];
			iss >> improperdihedralbond.global_ids[0] >> improperdihedralbond.global_ids[1] >> improperdihedralbond.global_ids[2] >> improperdihedralbond.global_ids[3] >> atomtypes[0] >> atomtypes[1] >> atomtypes[2] >> atomtypes[3] >> improperdihedralbond.psi0 >> improperdihedralbond.kpsi;
			improperdihedralbonds.entries.emplace_back(improperdihedralbond);
			break;
		}
		}
	}
}

SimulationFilesCollection::SimulationFilesCollection(const fs::path& workDir) {
	grofile = std::make_unique<ParsedGroFile>(workDir / "molecule/conf.gro");
	topfile = std::make_unique<ParsedTopologyFile>(workDir / "molecule/topol.top");
}
