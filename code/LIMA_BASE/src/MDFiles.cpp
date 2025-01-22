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
	auto constexpr trimSpaces = [](const std::string& str) -> std::string {
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
