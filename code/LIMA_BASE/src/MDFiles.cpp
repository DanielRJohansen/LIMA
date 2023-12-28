#include "MDFiles.h"



using namespace Filehandler;

GroRecord parseGroLine(const std::string& line) {
	GroRecord record;

	// Parse residue number (5 positions, integer)
	std::istringstream(line.substr(0, 5)) >> record.residue_number;

	// Parse residue name (5 characters)
	record.residue_name = line.substr(5, 5);
	record.residue_name.erase(std::remove(record.residue_name.begin(), record.residue_name.end(), ' '), record.residue_name.end());

	// Parse atom name (5 characters)
	record.atom_name = line.substr(10, 5);
	record.atom_name.erase(std::remove(record.atom_name.begin(), record.atom_name.end(), ' '), record.atom_name.end());

	// Parse atom number (5 positions, integer)
	std::istringstream(line.substr(15, 5)) >> record.gro_id;

	// Parse position (in nm, x y z in 3 columns, each 8 positions with 3 decimal places)
	std::istringstream(line.substr(20, 8)) >> record.position.x;
	std::istringstream(line.substr(28, 8)) >> record.position.y;
	std::istringstream(line.substr(36, 8)) >> record.position.z;

	if (line.size() >= 68) {
		// Parse velocity (in nm/ps (or km/s), x y z in 3 columns, each 8 positions with 4 decimal places)
		std::istringstream(line.substr(44, 8)) >> record.velocity.x;
		std::istringstream(line.substr(52, 8)) >> record.velocity.y;
		std::istringstream(line.substr(60, 8)) >> record.velocity.z;
	}

	return record;
}


ParsedGroFile MDFiles::loadGroFile(const std::filesystem::path& path) {
	SimpleParsedFile parsedfile = parseGroFile(path.string(), false);	// This step is slow because we read a stringstream and create strings, we could make the variables directly
	
	ParsedGroFile grofile{};

	for (const auto& row : parsedfile.rows) {
		if (row.section == "atoms") {
			grofile.atoms.push_back(parseGroLine(row.words[0]));
		}
		else if (row.section == "title") {
			grofile.title = row.words[0];
		}
		else if (row.section == "n_atoms") {
			grofile.n_atoms = stoi(row.words[0]);
		}
		else if (row.section == "box_size") {
			grofile.box_size = Float3{ stof(row.words[0]), stof(row.words[1]) ,stof(row.words[2]) };
		}
	}

	return grofile;
}

enum TopologySection {title, molecules, moleculetype, atoms, bonds, pairs, angles, dihedrals, impropers, no_section};

TopologySection getNextTopologySection(TopologySection current_section, const std::string& directive) {
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
		if (current_section != dihedrals)
			return dihedrals;
		else
			return impropers;
	}
	else {
		throw std::runtime_error(std::format("Got unexpected topology directive: {}", directive));
	}
}

std::string extractSectionName(const std::string& line) {
	size_t start = line.find('[');
	size_t end = line.find(']', start);
	if (start != std::string::npos && end != std::string::npos) {
		// Extract the text between '[' and ']'
		std::string sectionName = line.substr(start + 1, end - start - 1);
		return sectionName;
	}
	return "";
}

ParsedTopologyFile MDFiles::loadTopologyFile(const std::filesystem::path& path) {
	assert(path.extension().string() == std::string{ ".top" } || path.extension().string() == ".itp");

	ParsedTopologyFile topfile{};

	std::ifstream file;
	file.open(path);
	if (!file.is_open() || file.fail()) {
		throw std::runtime_error(std::format("Failed to open file {}\n", path.string()));
	}
	
	TopologySection current_section{ title };

	// Forward declaring for optimization reasons
	std::string line{}, word{};
	while (getline(file, line)) {

		//replaceTabs(line);

		std::vector<char> ignores = { ';', '#' };
		char delimiter = ' ';

		if (line.size() == 0) {
			current_section = no_section;
			continue;
		}

		if (line[0] == '[')
			current_section = getNextTopologySection(current_section, extractSectionName(line));

		if (line[0] == ';' && current_section != title) { continue; }	// Only title-sections reads the comments


		std::istringstream iss(line);
		switch (current_section)
		{
		case title:
			topfile.title.append(line);
			break;
		case molecules:
			throw std::runtime_error("This is not yet implemented");
			break;
		case moleculetype:
		{
			ParsedTopologyFile::MoleculetypeEntry moleculetype{};
			iss >> moleculetype.name >> moleculetype.nrexcl;
			topfile.moleculetypes.entries.emplace_back(moleculetype);
			break;
		}
		case atoms:
		{
			ParsedTopologyFile::AtomsEntry atom;
			iss >> atom.nr >> atom.type >> atom.resnr >> atom.residue >> atom.atom >> atom.cgnr >> atom.charge >> atom.mass;
			topfile.atoms.entries.emplace_back(atom);
			break;
		}
		case bonds: {
			ParsedTopologyFile::SingleBond singlebond{};
			iss >> singlebond.atom_indexes[0] >> singlebond.atom_indexes[1] >> singlebond.funct;
			topfile.singlebonds.entries.emplace_back(singlebond);
			break;
		}
		case pairs: {
			ParsedTopologyFile::Pair pair{};
			iss >> pair.atom_indexes[0] >> pair.atom_indexes[1] >> pair.funct;
			topfile.pairs.entries.emplace_back(pair);
			break;
		}
		case angles: {
			ParsedTopologyFile::AngleBond angle{};
			iss >> angle.atom_indexes[0] >> angle.atom_indexes[1] >> angle.atom_indexes[2] >> angle.funct;
			topfile.anglebonds.entries.emplace_back(angle);
			break; 
		}
		case dihedrals: {
			ParsedTopologyFile::DihedralBond dihedral{};
			iss >> dihedral.atom_indexes[0] >> dihedral.atom_indexes[1] >> dihedral.atom_indexes[2] >> dihedral.atom_indexes[3] >> dihedral.funct;
			topfile.dihedralbonds.entries.emplace_back(dihedral);
			break;
		}
		case impropers: {
			ParsedTopologyFile::ImproperDihedralBond improper{};
			iss >> improper.atom_indexes[0] >> improper.atom_indexes[1] >> improper.atom_indexes[2] >> improper.atom_indexes[3] >> improper.funct;
			topfile.improperdihedralbonds.entries.emplace_back(improper);
			break;
		}
		default:
			throw std::runtime_error("Illegal state");
			break;
		}
	}
	return topfile;
}