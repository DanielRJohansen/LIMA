#include "MDFiles.h"
#include "Format.h"
#include "Filehandling.h"
#include "MDFilesSerialization.h"
#include "ParallelFor.h"

#include <algorithm>
#include <charconv>
#include <cstdio>
#include <format>


using namespace FileUtils;
using namespace MDFiles;
namespace fs = std::filesystem;

namespace {
	std::string_view TrimSpaces(std::string_view sv) {
		while (!sv.empty() && (sv.front() == ' ' || sv.front() == '\t')) sv.remove_prefix(1);
		while (!sv.empty() && (sv.back() == ' ' || sv.back() == '\t')) sv.remove_suffix(1);
		return sv;
	}

	template <typename T>
	T ParseGroField(std::string_view line, size_t pos, size_t len) {
		const std::string_view field = TrimSpaces(line.substr(pos, len));
		T value{};
		const auto result = std::from_chars(field.data(), field.data() + field.size(), value);
		if (field.empty() || result.ec != std::errc{})
			throw std::runtime_error(Lima::Format("Failed to parse .gro line: \"{}\"", line));
		return value;
	}

	// Fixed width columns: resnr(5) resname(5) atomname(5) atomnr(5) pos(3x8) [vel(3x8)]
	GroRecord ParseGroLine(std::string_view line) {
		constexpr size_t minChars = 5 + 5 + 5 + 5 + 8 + 8 + 8;
		if (line.size() < minChars)
			throw std::runtime_error(Lima::Format("Too short .gro line: \"{}\"", line));

		GroRecord record;
		record.residue_number = ParseGroField<int>(line, 0, 5);
		record.residueName = TrimSpaces(line.substr(5, 5));
		record.atomName = TrimSpaces(line.substr(10, 5));
		record.gro_id = ParseGroField<int>(line, 15, 5);

		// Parse as double and then convert, same as the strtod we used to use, so the values are bit identical
		record.position = Float3{ ParseGroField<double>(line, 20, 8), ParseGroField<double>(line, 28, 8), ParseGroField<double>(line, 36, 8) };
		if (line.size() >= 68)	// [nm/ps]
			record.velocity = Float3{ ParseGroField<double>(line, 44, 8), ParseGroField<double>(line, 52, 8), ParseGroField<double>(line, 60, 8) };

		return record;
	}
}

// Fixed width columns, as ParseGroLine reads them. snprintf instead of an ostringstream per line, which made
// writing a 200k-atom file take over a second
static void AppendGroLine(std::string& out, const GroRecord& record) {
	const std::string_view residueName = record.residueName.View();
	const std::string_view atomName = record.atomName.View();
	char line[160];
	int length = std::snprintf(line, sizeof(line), "%-5d%-5.*s%-5.*s%5d%8.3f%8.3f%8.3f",
		record.residue_number, static_cast<int>(residueName.size()), residueName.data(),
		static_cast<int>(atomName.size()), atomName.data(), record.gro_id,
		record.position.x, record.position.y, record.position.z);
	if (record.velocity)
		length += std::snprintf(line + length, sizeof(line) - length, "%8.4f%8.4f%8.4f",
			record.velocity->x, record.velocity->y, record.velocity->z);
	out.append(line, length);
	out.push_back('\n');
}


GroFile::GroFile(const fs::path& path) : m_path(path){
	if (!(path.extension().string() == std::string{ ".gro" }))
		throw std::runtime_error("Expected .gro extension");
	if (!fs::exists(path))
		throw std::runtime_error(Lima::Format("File \"{}\" was not found", path.string()));

	lastModificationTimestamp = TimeSinceEpoch(fs::last_write_time(path));

	if (UseCachedBinaryFile(path)) {
		readGroFileFromBinaryCache(path, *this);
	}
	else {
		const std::string contents = ReadFileToString(path);

		// Split into lines, dropping the trailing \r of CRLF files
		std::vector<std::string_view> lines;
		lines.reserve(contents.size() / 45 + 3);	// Atom lines without velocities are 44 chars + newline
		for (std::string_view text = contents; !text.empty();) {
			const size_t newline = text.find('\n');
			std::string_view line = text.substr(0, newline);
			text.remove_prefix(newline == std::string_view::npos ? text.size() : newline + 1);
			if (!line.empty() && line.back() == '\r')
				line.remove_suffix(1);
			lines.push_back(line);
		}

		// Line 1 is the title, line 2 the atom count, then 1 line per atom, and finally the box
		if (lines.size() < 3)
			throw std::runtime_error(Lima::Format("File {} is too short to be a .gro file", path.string()));
		title = lines[0];

		size_t nAtoms = 0;
		const std::string_view countLine = TrimSpaces(lines[1]);
		if (std::from_chars(countLine.data(), countLine.data() + countLine.size(), nAtoms).ec != std::errc{})
			throw std::runtime_error(Lima::Format("Failed to read atom count in .gro file {}", path.string()));
		if (lines.size() < nAtoms + 3)
			throw std::runtime_error(Lima::Format(".gro file {} specifies {} atoms, but only has {} lines", path.string(), nAtoms, lines.size()));

		atoms.resize(nAtoms);
		ParallelUtils::ParallelForBlocked(nAtoms, [&](size_t i) { atoms[i] = ParseGroLine(lines[i + 2]); });

		// The box line has 3 or 9 values: the diagonal, then the off-diagonal ones of a triclinic box
		std::string_view boxLine = lines[nAtoms + 2];
		for (int dim = 0; dim < 3; dim++) {
			boxLine = TrimSpaces(boxLine);
			const auto result = std::from_chars(boxLine.data(), boxLine.data() + boxLine.size(), box_size[dim]);
			if (result.ec != std::errc{})
				throw std::runtime_error(Lima::Format("Failed to read box size in .gro file {}", path.string()));
			boxLine.remove_prefix(result.ptr - boxLine.data());
		}
		for (int offDiagonal = 0; offDiagonal < 6; offDiagonal++) {
			boxLine = TrimSpaces(boxLine);
			float value = 0.f;
			const auto result = std::from_chars(boxLine.data(), boxLine.data() + boxLine.size(), value);
			if (result.ec != std::errc{})
				break;
			triclinic |= value != 0.f;
			boxLine.remove_prefix(result.ptr - boxLine.data());
		}

		// Save a binary cached version of the file to we can read it faster next time
		WriteFileToBinaryCache(*this);
	}
}

void GroFile::printToFile(const std::filesystem::path& path) const {
	if (path.extension().string() != ".gro") { throw std::runtime_error(Lima::Format("Got {} extension, expected .gro", path.extension().string())); }
	if (!path.parent_path().empty())
		fs::create_directories(path.parent_path());

	std::ofstream file(path);
	if (!file.is_open()) {
		throw std::runtime_error(Lima::Format("Failed to open file {}", path.string()));
	}

	// Print the title and number of atoms
	file << title << "\n";
	file << atoms.size() << "\n";

	std::string lines;
	lines.reserve(atoms.size() * 69);
	for (const auto& atom : atoms)
		AppendGroLine(lines, atom);
	file.write(lines.data(), static_cast<std::streamsize>(lines.size()));

	// Print the box size
	file << box_size.x << " " << box_size.y << " " << box_size.z << "\n";

	file.close();

	// Also cache the file
	WriteFileToBinaryCache(*this, path);
}











// TODO: Remove this
SimulationFilesCollection::SimulationFilesCollection(const fs::path& workDir) {
	grofile = std::make_unique<GroFile>(workDir / "molecule/conf.gro");
	topfile = std::make_unique<TopologyFile>(workDir / "molecule/topol.top");
}


PDBfile::PDBfile(const fs::path& path) : mPath(path) {
	if (!(path.extension().string() == std::string{ ".pdb" }))
		throw std::runtime_error("Expected .pdb extension");
	if (!fs::exists(path))
		throw std::runtime_error(Lima::Format("File \"{}\" was not found", path.string()));



	std::ifstream file;
	file.open(path);
	if (!file.is_open() || file.fail()) {
		throw std::runtime_error(Lima::Format("Failed to open file {}\n", path.string()).c_str());
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
				ATOMS.emplace_back(atom);
			}
			else if (isHETATM) {
				HETATMS.emplace_back(atom);
			}
		}
	}
	
}
