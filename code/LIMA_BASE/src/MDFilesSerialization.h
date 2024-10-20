//#include "pch.h"

#include "MDFiles.h"
#include "Filehandling.h"

#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/optional.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/cereal.hpp>

#include <format>
#include <fstream>

using namespace FileUtils;
using namespace MDFiles;
namespace fs = std::filesystem;

inline int64_t TimeSinceEpoch(std::filesystem::file_time_type fileTime) {
	// Convert file_time_type to system_clock time_point
	auto sysTime = std::chrono::time_point_cast<std::chrono::system_clock::duration>(
		fileTime - std::filesystem::file_time_type::clock::now() + std::chrono::system_clock::now()
	);

	// Convert to milliseconds since epoch
	return std::chrono::duration_cast<std::chrono::milliseconds>(
		sysTime.time_since_epoch()
	).count();
}

inline constexpr uint64_t CacheVersionNumberValue() {
	const int cacheVersionNumber = 21;	// Modify this value each time we want to invalidate cached files made by previous versions of the program. 
	return 0xF0F0F0F0'00000000 + cacheVersionNumber;	// Cant just have a bunch of zeroes preceding the version, then we can't tell if the file is corrupted or not
}  


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
		archive(f.residueName);
		archive(f.atomName);
		archive(f.gro_id);
		archive(f.position);
		archive(f.velocity);
	}

	

	template <class Archive>
	void serialize(Archive& archive, TopologyFile::AtomsEntry& atom) {
		archive(atom.section_name);
		archive(atom.id);
		archive(atom.type);
		archive(atom.resnr);
		archive(atom.residue);
		archive(atom.atomname);
		archive(atom.cgnr);
		archive(atom.charge);
		archive(atom.mass);
	}
	template <class Archive, size_t N>
	void serialize(Archive& archive, TopologyFile::GenericBond<N>& bond) {
		archive(bond.ids);
		archive(bond.funct);
	}

	//template <class Archive>
	//void serialize(Archive& archive, TopologyFile::MoleculeEntry& molecule) {
	//	archive(molecule.name);
	//}

	template <class Archive>
	void serialize(Archive& archive, TopologyFile::Moleculetype& moleculetype) {
		archive(moleculetype.name);
		archive(moleculetype.nrexcl);
	}

	//template <class Archive, typename EntryType>
	//void serialize(Archive& archive, TopologyFile::Section<EntryType>& section) {
	//	archive(section.entries);
	//}

	template <class Archive>
	void serialize(Archive& ar, fs::path& p) {
		std::string path_str = p.string();
		ar(path_str);  // Serialize as a string
		if (Archive::is_loading::value) {
			p = fs::path(path_str);  // Convert back to fs::path when loading
		}
	}

	template <class Archive>
	void serialize(Archive& archive, TopologyFile::ForcefieldInclude& include) {
		archive(include.name);
	}

	template <class Archive>
	void serialize(Archive& archive, std::array<std::string, 2>& strArr) {
		archive(strArr[0]);
		archive(strArr[1]);
	}
} // namespace cereal


inline void readGroFileFromBinaryCache(const fs::path& path, GroFile& file) {
	std::ifstream is(path.string() + ".bin", std::ios::binary);
	if (!is.is_open()) {
		throw std::runtime_error("Failed to open file for reading");
	}

	cereal::BinaryInputArchive archive(is);
	uint64_t versionNumberValue;
	archive(versionNumberValue);	// ignore
	archive(file.lastModificationTimestamp);

	archive(file.title);
	archive(file.atoms);
	archive(file.box_size);

	file.readFromCache = true;
}
inline void WriteFileToBinaryCache(const GroFile& file, std::optional<fs::path> _path = std::nullopt) {
	const fs::path path = _path.value_or(file.m_path);
	if (path.empty())
		throw std::runtime_error("Tried to cache a Gro file with no path");
	std::ofstream os(path.string() + ".bin", std::ios::binary);
	if (!os.is_open()) {
		throw std::runtime_error("Failed to open file for writing: " + path.string() + ".bin");
	}

	cereal::BinaryOutputArchive archive(os);
	archive(CacheVersionNumberValue());
	archive(file.lastModificationTimestamp);

	archive(file.title);
	archive(file.atoms);
	archive(file.box_size);
}

//inline void readTopFileFromBinaryCache(const fs::path& path, TopologyFile::Moleculetype& moleculetype) {
//	std::ifstream is(path.string() + ".bin", std::ios::binary);
//	if (!is.is_open()) {
//		throw std::runtime_error("Failed to open file for reading");
//	}
//
//	cereal::BinaryInputArchive archive(is);
//	uint64_t versionNumberValue;
//	archive(versionNumberValue);	// ignore
//	archive(moleculetype.lastModificationTimestamp);
//
//	archive(moleculetype.name);
//	archive(moleculetype.nrexcl);
//
//	archive(moleculetype.atoms);
//	archive(moleculetype.singlebonds);
//	archive(moleculetype.pairs);
//	archive(moleculetype.anglebonds);
//	archive(moleculetype.dihedralbonds);
//	archive(moleculetype.improperdihedralbonds);
//
//	moleculetype.readFromCache = true;
//}
//
//inline void WriteFileToBinaryCache(const TopologyFile::Moleculetype& moleculetype, fs::path& dir) {
//	const fs::path path = dir / (moleculetype.name + ".top.bin");
//	if (dir.empty() || moleculetype.name.empty())
//		throw std::runtime_error("Tried to cache a Top file with no path");
//	{
//		std::ofstream os(path.string() + ".bin", std::ios::binary);
//		if (!os.is_open()) {
//			throw std::runtime_error("Failed to open file for writing: " + path.string() + ".bin");
//		}
//
//		cereal::BinaryOutputArchive archive(os);
//		archive(CacheVersionNumberValue());
//		archive(moleculetype.lastModificationTimestamp);
//
//		archive(moleculetype.name);
//		archive(moleculetype.nrexcl);
//
//		archive(moleculetype.atoms);
//		archive(moleculetype.singlebonds);
//		archive(moleculetype.pairs);
//		archive(moleculetype.anglebonds);
//		archive(moleculetype.dihedralbonds);
//		archive(moleculetype.improperdihedralbonds);
//	}
//}

// Checks if we have a valid cached binary of the file
static bool UseCachedBinaryFile(const fs::path& path) {
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
	if (versionNumber != CacheVersionNumberValue())
		return false;



	// Check if the file has been modified since the binary was created
	// If we are on linux and in the resources dir, we will NOT do the check, since the files were made on windows
	// and may mismatch AND we may not write to the resources dir
	// Actually lets just just not do the check at all in the resources dir, if we want to change these files
	// we can just delete the binaries, and they will be remade automatically
	if (path.parent_path().parent_path().parent_path() == FileUtils::GetLimaDir())
		return true;
#ifdef __linux__	
#endif	

	int64_t binaryModificationTime;
	file.read(reinterpret_cast<char*>(&binaryModificationTime), sizeof binaryModificationTime);
	const int64_t fileTime = TimeSinceEpoch(fs::last_write_time(path));

	// Check if the difference is less than 100 millisecond 
	return std::abs(fileTime - binaryModificationTime) < 100;
}
