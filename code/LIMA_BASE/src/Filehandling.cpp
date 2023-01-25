#include "Filehandling.h"



using namespace std;

bool Filehandler::ignoreRow(const vector<char>& ignores, const string& line) {
	if (line.length() == 0)
		return true;
	for (char c : ignores) {
		if (line[0] == c)
			return true;
	}
	return false;
}

bool Filehandler::ignoreWord(const vector<string>& ignores, const string& word) {
	if (word.length() == 0)
		return true;
	for (auto elem : ignores) {
		if (word == elem)
			return true;
	}
	return false;
}


// Uses ';' and ' ' as delimiters
vector<vector<string>> Filehandler::readFile(const string path, vector<char> comment_markers, std::vector<string> ignores, int end_at, bool verbose) {
	if (verbose) { cout << "Reading file " << path << endl; }
	fstream file;
	file.open(path);




	vector<vector<string>> rows;
	int row_cnt = 0;
	int ignore_cnt = 0;

	// Forward declaring for optimization reasons
	string line{}, element{}, element_nested{};
	while (getline(file, line)) {
		if (ignoreRow(comment_markers, line)) {
			ignore_cnt++;
			continue;
		}

		vector<string> row;
		stringstream ss(line);
		while (getline(ss, element, ' ')) {
			stringstream ss2 = stringstream(element);
			while (getline(ss2, element_nested, ';')) { 

				if (!ignoreWord(ignores, element_nested)) {
					row.push_back(element_nested);
				}
			}
		}

		rows.push_back(std::move(row));
		row_cnt++;

		if (row_cnt >= end_at)
			break;
	}

	if (verbose) {
		printf("%d rows read. %d rows ignored\n", row_cnt, ignore_cnt);
	}

	return rows;
}

map<string, double> Filehandler::parseINIFile(const string path) {
	// TODO: add read here
	//if (verbosity_level >= V1) { cout << "Reading particles from file " << path << "\n"; }
	fstream file;
	file.open(path);

	map<string, double> dict;

	string line;
	while (getline(file, line)) {
		int i = 0;
		string pair[2];
		stringstream ss(line);
		string word;
		while (getline(ss, word, '=')) {
			word.erase(std::remove_if(word.begin(), word.end(), std::isspace), word.end());
			pair[i++] = word;

			if (i == 2) { dict[pair[0]] = stod(pair[1]); }
		}
	}

	file.close();

	return dict;
}