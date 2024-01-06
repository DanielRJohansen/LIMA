// For generic file utilities, for specialized features use MDFiles.h instead
#pragma once



#include <iostream>
#include <vector>
#include <string>

#include <map>
#include <cassert>
#include <cstdint>
#include <limits>

struct SimpleParsedFile {
	struct Row {
		std::string section;
		std::vector<std::string> words;
	};

	std::vector<Row> rows;
};




// TODO: Why the fuck can i not make this a namespace???!
namespace Filehandler {
	static bool ignoreWord(const std::vector<std::string>& ignores, const std::string& word);

	static std::string pathJoin(std::string a, std::string b) { return a + "/" + b; }

	// Dunno if this works for folders too
	void assertPath(const std::string& path);

	bool fileExists(const std::string& path);

	void removeWhitespace(std::string& str);

	bool firstNonspaceCharIs(const std::string& str, char query);

	//void replaceTabs(std::string& str);

	std::string extractFilename(const std::string& path);

	std::map<std::string, double> parseINIFile(const std::string path);

	static std::vector<std::vector<std::string>> readFile(const std::string path, 
		std::vector<char> comment_markers = { ';', '/' },
		std::vector<std::string> ignore_words = { " "},
		int end_at = std::numeric_limits<int>::max(), bool verbose = false);

	static SimpleParsedFile parseItpFile(const std::string& path, bool verbose=true);
	SimpleParsedFile parseTopFile(const std::string& path, bool verbose);
	SimpleParsedFile parseLffFile(const std::string& path, bool verbose);
	SimpleParsedFile parsePrmFile(const std::string& path, bool verbose);
	SimpleParsedFile parseGroFile(const std::string& path, bool verbose);


	void createDefaultSimFilesIfNotAvailable(const std::string& dir, float boxsize_nm);	// creates conf topol and sim_params


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
};


