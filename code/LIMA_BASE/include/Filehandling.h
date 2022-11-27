#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>

#include <climits>
#include <map>

//using namespace std;





class Filehandler {
public:
	static bool ignoreRow(const std::vector<char>& ignores, const std::string& line);
	static bool ignoreWord(const std::vector<std::string>& ignores, const std::string& word);

	static std::map<std::string, double> parseINIFile(const std::string path);

	static std::vector<std::vector<std::string>> readFile(const std::string path, 
		std::vector<char> comment_markers = { ';', '/' },
		std::vector<std::string> ignore_words = { " "},
		int end_at = INT_MAX, bool verbose = false);
};

