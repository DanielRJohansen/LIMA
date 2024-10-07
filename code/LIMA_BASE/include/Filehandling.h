// For generic file utilities, for specialized features use MDFiles.h instead
#pragma once

#include <iostream>
#include <vector>
#include <string>

#include <unordered_map>
#include <cassert>
#include <cstdint>
#include <limits>

#include <filesystem>
#include <fstream>

struct SimpleParsedFile {
	struct Row {
		std::string section;
		std::vector<std::string> words;
	};

	std::vector<Row> rows;
};





namespace Filehandler {
	namespace fs = std::filesystem;

	static bool ignoreWord(const std::vector<std::string>& ignores, const std::string& word);

	// Dunno if this works for folders too
	void assertPath(const std::string& path);

	bool fileExists(const std::string& path);

	void removeWhitespace(std::string& str);

	bool firstNonspaceCharIs(const std::string& str, char query);

	//void replaceTabs(std::string& str);

	std::string extractFilename(const std::string& path);

	std::unordered_map<std::string, std::string> parseINIFile(const std::string& path, bool forceLowercase=false);

	static std::vector<std::vector<std::string>> readFile(const std::string path, 
		std::vector<char> comment_markers = { ';', '/' },
		std::vector<std::string> ignore_words = { " "},
		int end_at = std::numeric_limits<int>::max(), bool verbose = false);

	static SimpleParsedFile parseItpFile(const std::string& path, bool verbose=true);
	SimpleParsedFile parsePrmFile(const std::string& path, bool verbose);


	void createDefaultSimFilesIfNotAvailable(const std::string& dir, float boxsize_nm);	// creates conf topol and sim_params

	// Return the top level LIMA dir
	fs::path GetLimaDir();

	// These should be in interface maybe?
	template <typename T>
	static void dumpToFile(T* data, uint64_t n_datapoints, std::string file_path_s) {
#ifndef __linux__

		char* file_path;
		file_path = &file_path_s[0];

		const std::string str = std::to_string((long double)sizeof(T) * n_datapoints * 1e-6);
		//m_logger.print("Writing " + str + "MB to binary file " + file_path + "\n");

		FILE* file;

		if (!fopen_s(&file, file_path, "wb")) {

			assert(sizeof(T));
			assert(n_datapoints);

			fwrite(data, sizeof(T), n_datapoints, file);
			fclose(file);
		}
#else
		//file = fopen(file_path, "wb");
#endif
	}

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
};


