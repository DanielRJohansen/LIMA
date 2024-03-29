#include "MDFiles.h"

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

std::string ParsedTopologyFile::MoleculetypeEntry::composeString() const {
	std::ostringstream oss;
	oss << std::right << std::setw(10) << name << std::setw(10) << nrexcl;
	return oss.str();
}

std::string ParsedTopologyFile::AtomsEntry::composeString() const {
	std::ostringstream oss;
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
	return oss.str();
}






ParsedGroFile::ParsedGroFile(const fs::path& path) : m_path(path){
	assert(path.extension() == ".gro");



	std::ifstream file;
	file.open(path);
	if (!file.is_open() || file.fail()) {
		throw std::runtime_error(std::format("Failed to open file {}\n", path.string()).c_str());
	}



	int skipCnt = 2;	// First 2 lines are title and atom count

	const int min_chars = 5 + 5 + 5 + 5 + 8 + 8 + 8;

	// Forward declaring for optimization reasons
	std::string line{}, word{};
	while (getline(file, line)) {

		if (skipCnt > 0) {
			if (skipCnt == 2) {
				// 1st line is title
				title = line;
			}
			if (skipCnt == 1) {
				// 2nd line is atom count
				n_atoms = std::stoi(line);
				atoms.reserve(n_atoms);
			}

			skipCnt--;
			continue;
		}

		if (line.length() >= min_chars) {
			atoms.emplace_back(parseGroLine(line));
		}
		else if (line.size() > 0) {
			// Second last line, with box size
			int dim = 0;
			std::stringstream ss(line);
			while (std::getline(ss, word, ' ')) {
				if (!word.empty()) {
					assert(dim < 3);
					*box_size.placeAt(dim++) = std::stof(word);
				}
			}
		}
	}

	assert(atoms.size() == n_atoms);
}

enum TopologySection {title, molecules, moleculetype, atoms, bonds, pairs, angles, dihedrals, impropers, position_restraints, _system, cmap, no_section};

class TopologySectionGetter {
	int dihedralCount = 0;

public:
	TopologySection operator()(const std::string& directive) {
		if (directive == "molecules") {
			return molecules;
		}
		else if (directive == "moleculetype") {
			return moleculetype;
		}
		else if (directive == "atoms") {
			return atoms;
		}
		else if (directive == "bonds") {
			return bonds;
		}
		else if (directive == "pairs") {
			return pairs;
		}
		else if (directive == "angles") {
			return angles;
		}
		else if (directive == "dihedrals") {
			dihedralCount++;
			if (dihedralCount == 1) {
				return dihedrals;
			}
			else if (dihedralCount == 2) {
				return impropers;
			}
			else
				throw std::runtime_error("Encountered the 'dihedral' directive more than 2 times in the same .itp/.top file");
		}
		else if (directive == "position_restraints") {
			return position_restraints;
		}
		else if (directive == "system") {
			return _system;
		}
		else if (directive == "cmap") {
			return cmap;
		}
		else {
			throw std::runtime_error(std::format("Got unexpected topology directive: {}", directive));
		}
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

bool containsSubString(const std::istringstream& stream, const std::string& query) {
	std::istringstream issCopy(stream.str()); // Create a copy of the original stream
	std::string word;
	while (issCopy >> word) {
		if (word == query) {
			return true;
		}
	}
	return false;
}



std::unique_ptr<ParsedTopologyFile> MDFiles::loadTopologyFile(const fs::path& path) {
	assert(path.extension().string() == std::string{ ".top" } || path.extension().string() == ".itp");

	std::unique_ptr<ParsedTopologyFile> topfile = std::make_unique<ParsedTopologyFile>();

	topfile->m_path = path;

	const char delimiter = ' ';
	const std::vector<char> ignores = { ';', '#' };


	std::ifstream file;
	file.open(path);
	if (!file.is_open() || file.fail()) {
		throw std::runtime_error(std::format("Failed to open file {}\n", path.string()));
	}
	
	TopologySection current_section{ title };
	TopologySectionGetter getTopolSection{};

	// Forward declaring for optimization reasons
	std::string line{}, word{};

	std::string sectionname = "";


	//std::vector<ParsedTopologyFile> includeTopologies;

	std::vector<std::future<std::unique_ptr<ParsedTopologyFile>>> includeTopologies;

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
		if (firstNonspaceCharIs(line, ';') && current_section != title && current_section != atoms) { continue; }	// Only title-sections + atoms reads the comments


		std::istringstream iss(line);
		switch (current_section)
		{
		case title:
			topfile->title.append(line + "\n");	// +\n because getline implicitly strips it away.
			break;
		case molecules: {
				// TODO: This is incorrect, we should have a section where we look for "#include" and include those files
			std::string include_name;
			iss >> include_name;

			// Handle the case where there simply a random SOL that does not refer to a file.. Terrible standard...
			if (include_name == "SOL")
				continue;

			fs::path include_topol_path = path;
			include_topol_path.replace_filename("topol_" + include_name + ".itp");
			//ParsedTopologyFile* f = loadTopologyFile(include_topol_path).release();
			//topfile->molecules.entries.emplace_back(ParsedTopologyFile::MoleculeEntry{ include_name, f });

			includeTopologies.emplace_back(std::async(std::launch::async, &MDFiles::loadTopologyFile, include_topol_path));

			break;
		}
		case moleculetype:
		{
			ParsedTopologyFile::MoleculetypeEntry moleculetype{};
			iss >> moleculetype.name >> moleculetype.nrexcl;
			topfile->moleculetypes.entries.emplace_back(moleculetype);
			break;
		}
		case atoms:
		{
			if (firstNonspaceCharIs(line, ';') ) {	

				// TODO: Test for residue or lipid_section in the [1] position of the comment instead
				if (!containsSubString(iss, "type"))// We want any comment that is NOT the legend and the legend always contains "type"
					sectionname = line;
			}
			else {
				ParsedTopologyFile::AtomsEntry atom;
				iss >> atom.nr >> atom.type >> atom.resnr >> atom.residue >> atom.atomname >> atom.cgnr >> atom.charge >> atom.mass;
				topfile->atoms.entries.emplace_back(atom);

				if (sectionname != "") {
					topfile->atoms.entries.back().section_name = sectionname;
					sectionname = "";
				}
			}
			break;
		}
		case bonds: {
			ParsedTopologyFile::SingleBond singlebond{};
			iss >> singlebond.atom_indexes[0] >> singlebond.atom_indexes[1] >> singlebond.funct;
			topfile->singlebonds.entries.emplace_back(singlebond);
			break;
		}
		case pairs: {
			ParsedTopologyFile::Pair pair{};
			iss >> pair.atom_indexes[0] >> pair.atom_indexes[1] >> pair.funct;
			topfile->pairs.entries.emplace_back(pair);
			break;
		}
		case angles: {
			ParsedTopologyFile::AngleBond angle{};
			iss >> angle.atom_indexes[0] >> angle.atom_indexes[1] >> angle.atom_indexes[2] >> angle.funct;
			topfile->anglebonds.entries.emplace_back(angle);
			break; 
		}
		case dihedrals: {
			ParsedTopologyFile::DihedralBond dihedral{};
			iss >> dihedral.atom_indexes[0] >> dihedral.atom_indexes[1] >> dihedral.atom_indexes[2] >> dihedral.atom_indexes[3] >> dihedral.funct;
			topfile->dihedralbonds.entries.emplace_back(dihedral);
			break;
		}
		case impropers: {
			ParsedTopologyFile::ImproperDihedralBond improper{};
			iss >> improper.atom_indexes[0] >> improper.atom_indexes[1] >> improper.atom_indexes[2] >> improper.atom_indexes[3] >> improper.funct;
			topfile->improperdihedralbonds.entries.emplace_back(improper);
			break;
		}
		default:
			// Do nothing
			//throw std::runtime_error("Illegal state");
			break;
		}
	}


	// Add all include topologies from other threads
	for (auto& includeTopologyFuture: includeTopologies) {
		try {
			std::unique_ptr<ParsedTopologyFile> includeTopology = includeTopologyFuture.get(); // This will block until the future is ready

			if (includeTopology) {
				topfile->molecules.entries.emplace_back(ParsedTopologyFile::MoleculeEntry{ includeTopology->m_path.string(), includeTopology.release()});
			}
		}
		catch (const std::exception& e) {
			std::cerr << "Failed to load included topology file\n\t " << e.what() << std::endl;
			// Handle the error or fail according to your application's needs
		}
	}



	return topfile;
}

void ParsedGroFile::printToFile(const std::filesystem::path& path) const {
	if (path.extension().string() != ".gro") { throw std::runtime_error(std::format("Got {} extension, expected .gro", path.extension().string())); }


	std::ofstream file(path);
	if (!file.is_open()) {
		throw std::runtime_error(std::format("Failed to open file {}", path.string()));
	}

	// Print the title and number of atoms
	file << title << "\n";
	file << n_atoms << "\n";

	// Iterate over atoms and print them
	for (const auto& atom : atoms) {
		// You need to define how GroRecord is formatted
		file << composeGroLine(atom) << "\n";
	}

	// Print the box size
	file << box_size.x << " " << box_size.y << " " << box_size.z << "\n";

	file.close();
}

void ParsedGroFile::addEntry(std::string residue_name, std::string atom_name, const Float3& position) {
	if (n_atoms == 0) {
		atoms.push_back({ 0, residue_name, atom_name, 0, position, std::nullopt });
	}
	else {
		atoms.push_back({atoms.back().residue_number+1, residue_name, atom_name, atoms.back().gro_id+1, position, std::nullopt});
	}
	n_atoms++;
}

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

	//if (!molecules.entries.empty()) { file << molecules.composeString(); }
	if (!moleculetypes.entries.empty()) { file << moleculetypes.composeString(); }
	
	
	file << atoms.composeString();
	file << singlebonds.composeString();
	file << pairs.composeString();
	file << anglebonds.composeString();
	file << dihedralbonds.composeString();
	file << improperdihedralbonds.composeString();
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
	topfile = loadTopologyFile(workDir / "molecule/topol.top");
}