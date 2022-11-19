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

	static std::map<std::string, double> parseINIFile(std::string path);

	static std::vector<std::vector<std::string>> readFile(std::string path, int end_at = INT_MAX, bool verbose = false);
};

