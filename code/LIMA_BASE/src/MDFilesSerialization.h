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

	//template <class Archive>
	//void save(Archive& archive, const Float3& f) {
	//	archive(f.x);
	//	archive(f.y);
	//	archive(f.z);
	//}

	//template <class Archive>
	//void load(Archive& archive, Float3& f) {
	//	archive(f.x);
	//	archive(f.y);
	//	archive(f.z);
	//}

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


} // namespace cereal



void writeFileAsBinaryCache(const ParsedGroFile& file) {
	std::ofstream os(file.m_path.string() + ".bin", std::ios::binary);
	if (!os.is_open()) {
		throw std::runtime_error("Failed to open file for writing");
	}

	cereal::BinaryOutputArchive archive(os);
	archive(file.lastModificationTimestamp);
	archive(file.title);
	archive(file.atoms);
	archive(file.box_size);
}

ParsedGroFile readFileFromBinaryCache(const fs::path& path) {
	std::ifstream is(path.string() + ".bin", std::ios::binary);
	if (!is.is_open()) {
		throw std::runtime_error("Failed to open file for reading");
	}
	ParsedGroFile file;

	cereal::BinaryInputArchive archive(is);
	archive(file.lastModificationTimestamp);
	archive(file.title);
	archive(file.atoms);
	archive(file.box_size);	
	
	file.readFromCache = true;
	return std::move(file);
}

// Checks if we have a valid cached binary of the file
static bool UseCachedBinaryFile(const fs::path& path, fs::file_time_type requredModificationTime) {
	const fs::path binaryPath = path.string() + ".bin";
	if (!fs::exists(binaryPath))
		return false;

	std::ifstream file;
	file.open(binaryPath);
	if (!file.is_open() || file.fail()) {
		throw std::runtime_error(std::format("Failed to open file {}\n", path.string()).c_str());
	}

	fs::file_time_type binaryModificationTime;;
	file.read(reinterpret_cast<char*>(&binaryModificationTime), sizeof(fs::file_time_type));

	return binaryModificationTime == requredModificationTime;
}
