#include "Filehandling.h"



using namespace std;

bool Filehandler::ignoreRow(const vector<char>& ignores, const string& line) {
	if (line.length() == 0)
		return true;
	for (auto& c : ignores) {
		if (line[0] == c)
			return true;
	}
	return false;
}

bool Filehandler::ignoreWord(const vector<string>& ignores, const string& word) {
	if (word.length() == 0)
		return true;
	for (auto& elem : ignores) {
		if (word == elem)
			return true;
	}
	return false;
}


// Uses ';' and ' ' as delimiters
vector<vector<string>> Filehandler::readFile(const string path, vector<char> comment_markers, std::vector<string> ignores, int end_at, bool verbose) {
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
		if (row.empty()) { continue; }	// This case happens when a line contains 1 or more spaces, but no words. Space are not regarded as comments, since the separate entries in a line

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

SimpleParsedFile Filehandler::parseItpFile(const std::string& path, bool verbose)
{
	fstream file;
	file.open(path);

	SimpleParsedFile parsedfile;

	const auto ignores = { ';', '#'};	// comment marker
	const auto delimiter = ' ';		// separator


	string current_section = "none";

	int row_cnt = 0;
	int ignore_cnt = 0;

	// Forward declaring for optimization reasons
	string line{}, word{};
	while (getline(file, line)) {

		// Check if line is a comment
		if (ignoreRow(ignores, line)) {
			ignore_cnt++;
			continue;
		}

		vector<string> row;
		stringstream ss(line);
		while (getline(ss, word, delimiter)) {
			if (!word.empty()) {
				row.push_back(word);
			}
			
		}

		if (row.empty()) { continue; }	// This case happens when a line contains 1 or more spaces, but no words. Space are not regarded as comments, since the separate entries in a line

		// Check for new section - expecting syntax like this:"[ angletypes ]"
		if (row.size() == 3 && row[0][0] == '[') {

			// Need to handle a special case, because some fuckwits used the same keyword twice - straight to jail!

			if (current_section == "dihedraltypes" && row[1] == "dihedraltypes") {	// Workaround for itp files
				current_section = "improperdihedraltypes";
			}
			else if (current_section == "dihedrals" && row[1] == "dihedrals") {		// Workaround for top files
				current_section = "improperdihedrals";
			}
			else {
				current_section = row[1];
			}
			
			continue;
		}

		parsedfile.rows.push_back({ current_section, row });

		row_cnt++;

	}

	if (verbose) {
		printf("%d rows read. %d rows ignored\n", row_cnt, ignore_cnt);
	}	
	
	return parsedfile;
}

SimpleParsedFile Filehandler::parseTopFile(const std::string& path, bool verbose)
{
	// TODO: make both of these functions call a generic parsefile with their own ignores, specialcase handling of sections and more
	return parseItpFile(path, verbose);
}
