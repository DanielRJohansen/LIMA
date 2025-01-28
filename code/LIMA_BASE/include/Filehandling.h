// For generic file utilities, for specialized features use MDFiles.h instead
#pragma once

#include "LimaTypes.cuh"

#include <cassert>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>


namespace FileUtils {
	namespace fs = std::filesystem;

	static bool ignoreWord(const std::vector<std::string>& ignores, const std::string& word);

	void removeWhitespace(std::string& str);

	bool firstNonspaceCharIs(const std::string& str, char query);

	std::unordered_map<std::string, std::string> parseINIFile(const std::string& path, bool forceLowercase=false);

	std::string_view ExtractBetweenQuotemarks(const std::string& input);

	// Return the top level LIMA dir
	fs::path GetLimaDir();

	std::vector<std::array<fs::path, 2>> GetAllGroItpFilepairsInDir(const fs::path& dir);

	template <typename T>
	std::vector<T> ReadBinaryFileIntoVector(const fs::path& path) {
		if (!fs::exists(path)) {
			throw std::runtime_error("File does not exist: " + path.string());
		}

		std::ifstream file(path, std::ios::binary | std::ios::ate);
			if (!file) {
				throw std::runtime_error("Failed to open file: " + path.string());
			}

		std::streamsize size = file.tellg();
			if (size % sizeof(T) != 0) {
				throw std::runtime_error("File size is not a multiple of type size.");
			}

		file.seekg(0, std::ios::beg);
		std::vector<T> buffer(size / sizeof(T));

		if (!file.read(reinterpret_cast<char*>(buffer.data()), size)) {
			throw std::runtime_error("Error reading file: " + path.string());
		}

		return buffer;
	}

	template <typename T>
	void WriteVectorToBinaryFile(const fs::path& path, const std::vector<T>& vec) {
		if (path.empty())
			throw std::runtime_error("Empty path");

		std::ofstream file(path, std::ios::binary);
		if (!file) {
			throw std::runtime_error("Failed to open file: " + path.string());
		}

		file.write(reinterpret_cast<const char*>(vec.data()), vec.size() * sizeof(T));
		if (!file) {
			throw std::runtime_error("Error writing to file: " + path.string());
		}
	}

	//Steps through the file untill it finds a #endif, throws if it reaches EoF
	void SkipIfdefBlock(std::ifstream& file);

	// Returns true if the caller should move on to the next line
	bool ChecklineForIfdefAndSkipIfFound(std::ifstream& file, const std::string& line, const std::unordered_set<std::string>& defines);

	std::optional<std::string> ChechlineForDefine(const std::string& line);

	std::vector<Float3> ReadCsvAsVectorOfFloat3(const fs::path& path);
};


