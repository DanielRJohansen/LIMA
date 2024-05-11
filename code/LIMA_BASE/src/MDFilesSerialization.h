//#include "pch.h"

#include "MDFiles.h"
#include "Filehandling.h"

#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/optional.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/cereal.hpp>

#include <fstream>

using namespace Filehandler;
using namespace MDFiles;
namespace fs = std::filesystem;


const int __CacheVersionNumber = 2;	// Modify this value each time we want to invalidate cached files made by previous versions of the program. 
const uint64_t CacheVersionNumberValue = 0xF0F0F0F0'00000000 + __CacheVersionNumber;	// Cant just have a bunch of zeroes preceding the version, then we can't tell if the file is corrupted or not


namespace cereal {

	template <class Archive>
	void save(Archive& archive, const std::filesystem::file_time_type& f) {
		char buffer[sizeof(std::filesystem::file_time_type)];
		static_assert(sizeof(buffer) == 8, "Expected file_time_type to be 8 bytes");
		std::memcpy(buffer, &f, sizeof(buffer));
		archive(cereal::binary_data(buffer, sizeof(buffer)));
	}

	template <class Archive>
	void load(Archive& archive, std::filesystem::file_time_type& f) {
		char buffer[sizeof(std::filesystem::file_time_type)];
		archive(cereal::binary_data(buffer, sizeof(buffer)));
		std::memcpy(&f, buffer, sizeof(buffer));
	}

	template <class Archive>
	void serialize(Archive& archive, Float3& f) {
		archive(f.x);
		archive(f.y);
		archive(f.z);
	}

	template <class Archive>
	void serialize(Archive& archive, GroRecord& f) {
		archive(f.residue_number);
		archive(f.residue_name);
		archive(f.atom_name);
		archive(f.gro_id);
		archive(f.position);
		archive(f.velocity);
	}

	

	template <class Archive>
	void serialize(Archive& archive, ParsedTopologyFile::AtomsEntry& atom) {
		archive(atom.section_name);
		archive(atom.nr);
		archive(atom.type);
		archive(atom.resnr);
		archive(atom.residue);
		archive(atom.atomname);
		archive(atom.cgnr);
		archive(atom.charge);
		archive(atom.mass);
	}
	template <class Archive, size_t N>
	void serialize(Archive& archive, ParsedTopologyFile::GenericBond<N>& bond) {
		archive(bond.atomGroIds);
		archive(bond.funct);
	}

	template <class Archive>
	void serialize(Archive& archive, ParsedTopologyFile::MoleculeEntry& molecule) {
		archive(molecule.name);
	}

	template <class Archive>
	void serialize(Archive& archive, ParsedTopologyFile::MoleculetypeEntry& moleculetype) {
		archive(moleculetype.name);
		archive(moleculetype.nrexcl);
	}

	template <class Archive, typename EntryType>
	void serialize(Archive& archive, Section<EntryType>& section) {
		archive(section.entries);
	}

} // namespace cereal


ParsedGroFile readGroFileFromBinaryCache(const fs::path& path) {
	std::ifstream is(path.string() + ".bin", std::ios::binary);
	if (!is.is_open()) {
		throw std::runtime_error("Failed to open file for reading");
	}
	ParsedGroFile file;

	cereal::BinaryInputArchive archive(is);
	uint64_t versionNumberValue;
	archive(versionNumberValue);	// ignore
	archive(file.lastModificationTimestamp);

	archive(file.title);
	archive(file.atoms);
	archive(file.box_size);

	file.readFromCache = true;
	return std::move(file);
}
void WriteFileToBinaryCache(const ParsedGroFile& file, std::optional<fs::path> _path = std::nullopt) {
	const fs::path path = _path.value_or(file.m_path);
	if (path.empty())
		throw std::runtime_error("Tried to cache a Gro file with no path");
	std::ofstream os(path.string() + ".bin", std::ios::binary);
	if (!os.is_open()) {
		throw std::runtime_error("Failed to open file for writing");
	}

	cereal::BinaryOutputArchive archive(os);
	archive(CacheVersionNumberValue);
	archive(file.lastModificationTimestamp);

	archive(file.title);
	archive(file.atoms);
	archive(file.box_size);
}

void readTopFileFromBinaryCache(const fs::path& path, ParsedTopologyFile& file) {
	std::ifstream is(path.string() + ".bin", std::ios::binary);
	if (!is.is_open()) {
		throw std::runtime_error("Failed to open file for reading");
	}

	cereal::BinaryInputArchive archive(is);
	uint64_t versionNumberValue;
	archive(versionNumberValue);	// ignore
	archive(file.lastModificationTimestamp);

	archive(file.title);
	archive(file.molecules);	// Only write the name of the file
	archive(file.moleculetypes);
	archive(file.GetLocalAtoms());
	archive(file.GetLocalSinglebonds());
	archive(file.GetLocalPairs());
	archive(file.GetLocalAnglebonds());
	archive(file.GetLocalDihedralbonds());
	archive(file.GetLocalImproperDihedralbonds());
	file.readFromCache = true;
}

void WriteFileToBinaryCache(const ParsedTopologyFile& file, std::optional<fs::path> _path = std::nullopt ) {
	const fs::path path = _path.value_or(file.path);
	if (path.empty())
		throw std::runtime_error("Tried to cache a Top file with no path");
	std::ofstream os(path.string() + ".bin", std::ios::binary);
	if (!os.is_open()) {
		throw std::runtime_error("Failed to open file for writing");
	}

	cereal::BinaryOutputArchive archive(os);
	archive(CacheVersionNumberValue);
	archive(file.lastModificationTimestamp);

	archive(file.title);
	archive(file.molecules);	// Only write the path of the file
	archive(file.moleculetypes);
	archive(file.GetLocalAtoms());
	archive(file.GetLocalSinglebonds());
	archive(file.GetLocalPairs());
	archive(file.GetLocalAnglebonds());
	archive(file.GetLocalDihedralbonds());
	archive(file.GetLocalImproperDihedralbonds());
}

// Checks if we have a valid cached binary of the file
static bool UseCachedBinaryFile(const fs::path& path, fs::file_time_type requredModificationTime) {
	if (!ENABLE_FILE_CACHING)
		return false;

	const fs::path binaryPath = path.string() + ".bin";
	if (!fs::exists(binaryPath))
		return false;


	// Check if the file exists and is at least 8 bytes
	if (fs::file_size(binaryPath) < (sizeof(fs::file_time_type) + sizeof(uint64_t)))
		return false;

	std::ifstream file;
	file.open(binaryPath);
	if (!file.is_open() || file.fail()) {
		throw std::runtime_error(std::format("Failed to open file {}\n", path.string()).c_str());
	}

	// Read the version number as a 64-bit unsigned integer
	uint64_t versionNumber;
	file.read(reinterpret_cast<char*>(&versionNumber), sizeof(uint64_t));
	if (versionNumber != CacheVersionNumberValue)
		return false;


	fs::file_time_type binaryModificationTime;;
	file.read(reinterpret_cast<char*>(&binaryModificationTime), sizeof(fs::file_time_type));
	return binaryModificationTime == requredModificationTime;
}