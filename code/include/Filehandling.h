#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>

#include <climits>

using namespace std;






class Filehandler {
public:
	static bool ignoreRow(vector<char> ignores, string line);

	static vector<vector<string>> readFile(string path, int end_at = INT_MAX, bool verbose = false);
};

