#pragma once

#include <iostream>
#include <vector>
#include <string>

#include <map>
#include <cassert>

struct SimpleParsedFile {
	struct Row {
		std::string section;
		std::vector<std::string> words;
	};

	std::vector<Row> rows;
};

// TODO: Why the fuck can i not make this a namespace???!
struct Filehandler {
	static bool ignoreWord(const std::vector<std::string>& ignores, const std::string& word);

	static std::string pathJoin(std::string a, std::string b) { return a + "/" + b; }

	// Dunno if this works for folders too
	static void assertPath(const std::string& path);

	static std::map<std::string, double> parseINIFile(const std::string path);

	static std::vector<std::vector<std::string>> readFile(const std::string path, 
		std::vector<char> comment_markers = { ';', '/' },
		std::vector<std::string> ignore_words = { " "},
		int end_at = INT_MAX, bool verbose = false);

	static SimpleParsedFile parseItpFile(const std::string& path, bool verbose=true);
	static SimpleParsedFile parseTopFile(const std::string& path, bool verbose);
	static SimpleParsedFile parseLffFile(const std::string& path, bool verbose);
	static SimpleParsedFile parsePrmFile(const std::string& path, bool verbose);

	// These should be in interface maybe?
	template <typename T>
	static void dumpToFile(T* data, uint64_t n_datapoints, std::string file_path_s) {
		char* file_path;
		file_path = &file_path_s[0];

		const std::string str = std::to_string((long double)sizeof(T) * n_datapoints * 1e-6);
		//m_logger.print("Writing " + str + "MB to binary file " + file_path + "\n");

		FILE* file;

#ifndef __linux__
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
};


