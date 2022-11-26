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
	static bool ignoreRow(std::vector<char> ignores, std::string line);
	static bool ignoreWord(std::vector<std::string> ignores, std::string word);

	static std::map<std::string, double> parseINIFile(std::string path);

	static std::vector<std::vector<std::string>> readFile(std::string path, 
		std::vector<char> comment_markers = { ';', '/' },
		std::vector<std::string> ignore_words = { " "},
		int end_at = INT_MAX, bool verbose = false);
};

