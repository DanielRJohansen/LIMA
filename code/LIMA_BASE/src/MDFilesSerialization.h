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
	void serialize(Archive& archive, ParsedTopologyFile::SingleBond& bond) {
		archive(bond.atomGroIds[0]);
		archive(bond.atomGroIds[1]);
		archive(bond.funct);
	}


	template <class Archive>
	void serialize(Archive& archive, ParsedTopologyFile::MoleculeEntry& molecule) {
		archive(molecule.name);
		// We can't serialize the other molecule recursively, so we only save the path, and then load from that path in the constructor
	}

	template <class Archive>
	void serialize(Archive& archive, ParsedTopologyFile::MoleculetypeEntry& moleculetype) {
		archive(moleculetype.name);
		archive(moleculetype.nrexcl);
		// We can't serialize the other molecule recursively, so we only save the path, and then load from that path in the constructor
	}

	template <class Archive, typename EntryType>
	void serialize(Archive& archive, Section<EntryType>& section) {
		archive(section.title);
		archive(section.entries);
	}

	//template <class Archive, typename EntryType>
	//void save(Archive& archive, const Section<EntryType>& section) {
	//	archive(section.title);
	//	archive(section.entries);
	//}



	//template <class Archive>
	//void load(Archive& archive, Section<ParsedTopologyFile::SingleBond>& section) {
	//	archive(section.title);
	//	archive(section.entries);
	//}


	//template <class Archive>
	//void serialize(Archive& archive, Section<ParsedTopologyFile::SingleBond>& section) {
	//	archive(section.title);
	//	//archive(section.entries);
	//}

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
void writeFileAsBinaryCache(const ParsedGroFile& file) {
	std::ofstream os(file.m_path.string() + ".bin", std::ios::binary);
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
	archive(file.molecules);	// Only write the path of the file
	archive(file.moleculetypes);
	archive(file.atoms);
	archive(file.singlebonds);
	archive(file.pairs);
	archive(file.anglebonds);
	archive(file.dihedralbonds);
	archive(file.improperdihedralbonds);
	file.readFromCache = true;
}

void WriteTopFileAsBinaryCache(const ParsedTopologyFile& file) {
	std::ofstream os(file.m_path.string() + ".bin", std::ios::binary);
	if (!os.is_open()) {
		throw std::runtime_error("Failed to open file for writing");
	}

	cereal::BinaryOutputArchive archive(os);
	archive(CacheVersionNumberValue);
	archive(file.lastModificationTimestamp);

	archive(file.title);
	archive(file.molecules);	// Only write the path of the file
	archive(file.moleculetypes);
	archive(file.atoms);
	archive(file.singlebonds);
	archive(file.pairs);
	archive(file.anglebonds);
	archive(file.dihedralbonds);
	archive(file.improperdihedralbonds);
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
