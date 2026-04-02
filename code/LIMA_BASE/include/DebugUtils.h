#pragma once

#include "Filehandling.h"





namespace DebugUtils {
	namespace fs = std::filesystem;
	template <typename T>
	bool VerifyIdentical(const std::vector<T>& data, const std::string& name) {
		//return true;
		const fs::path filePath = FileUtils::GetLimaDir() / "dev" / "etc" / (name + ".bin");

		if (fs::exists(filePath)) {
			const std::vector<T> fileData = FileUtils::ReadBinaryFileIntoVector<T>(filePath);

			if (data.size() != fileData.size()) {
				std::cout << name << " size mismatch for " << name << ": " << data.size() << " vs " << fileData.size() << std::endl;
				FileUtils::WriteVectorToBinaryFile(filePath, data);
				return false;
			}
			for (size_t i = 0; i < data.size(); ++i) {
				if (data[i] != fileData[i]) {
					std::cout << name << " data mismatch at index " << i << std::endl;
					FileUtils::WriteVectorToBinaryFile(filePath, data);
					return false;
				}
			}
		}

		FileUtils::WriteVectorToBinaryFile(filePath, data);
		return true;
	}

	template <typename T>
	bool VerifyIdentical(const T* const devPtr, size_t count, const std::string& name, int step) {
		if constexpr (DETERMINISTIC_CHECKS) {			
			return VerifyIdentical(GenericCopyToHost(devPtr, count), name + std::to_string(step));
		}
	}

}