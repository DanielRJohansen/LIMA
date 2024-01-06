#include "Filehandling.h"

#include <assert.h>
#include <algorithm>
#include <functional>
#include <map>
#include <array>
#include <fstream>
#include <filesystem>
#include <format>
#include <cctype>
#include <vector>
using std::string, std::vector, std::map, std::stringstream;


/// <summary>
/// Returns true for a new section. A new section mean the current line (title) is skipped
/// The function can also modify the skipCnt ref so the following x lines are skipped
/// </summary>
using SetSectionFunction = std::function<bool(const std::vector<string>& row, string& section, int& skipCnt)>;

bool ignoreRow(const vector<char>& ignores, const string& line) {
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

void Filehandler::assertPath(const std::string& path) {
	if (!std::filesystem::exists(path)) {
		std::cerr << "Could not find path: " << path << "\n";
		abort();
	}
}

bool Filehandler::fileExists(const std::string& path) {
	return std::filesystem::exists(path);
}



void Filehandler::removeWhitespace(std::string& str) {
	str.erase(std::remove_if(str.begin(), str.end(), 
		[](unsigned char ch) {return std::isspace(static_cast<int>(ch));}),
		str.end());
}

bool Filehandler::firstNonspaceCharIs(const std::string& str, char query) {
	auto first_non_space = std::find_if(str.begin(), str.end(), [](unsigned char ch) {
		return !std::isspace(static_cast<int>(ch));
	});

	return (first_non_space != str.end() && *first_non_space == query);
}

std::string Filehandler::extractFilename(const std::string& path) {
	// Find the last occurrence of both forward slash and backslash
	size_t lastSlashPos = path.find_last_of("/\\");

	// Check if a slash is found
	if (lastSlashPos != std::string::npos) {
		// Extract the substring from the last slash to the end of the string
		return path.substr(lastSlashPos + 1);
	}
	else {
		// If no slash is found, the entire path is the file name
		return path;
	}
}

// Uses ';' and ' ' as delimiters
vector<vector<string>> Filehandler::readFile(const string path, vector<char> comment_markers, std::vector<string> ignores, int end_at, bool verbose) {
	std::ifstream file;
	file.open(path);
	if (!file.is_open() || file.fail()) {
        throw std::runtime_error(std::format("Failed to open file {}\n", path).c_str());
    }



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
	std::ifstream file;
	file.open(path);
	if (!file.is_open() || file.fail()) {
        throw std::runtime_error(std::format("Failed to open file {}\n", path).c_str());
    }

	map<string, double> dict;

	string line;
	while (getline(file, line)) {
		int i = 0;
		string pair[2];
		stringstream ss(line);
		string word;
		while (getline(ss, word, '=')) {
			word.erase(std::remove_if(word.begin(), word.end(), isspace), word.end());
			pair[i++] = word;

			if (i == 2) { dict[pair[0]] = stod(pair[1]); }
		}
	}

	file.close();

	return dict;
}

void replaceTabs(std::string& str) {
	// Find the first occurrence of '\t' in the string
	size_t found = str.find('\t');

	// Continue replacing '\t' with spaces until no more occurrences are found
	while (found != std::string::npos) {
		// Replace the '\t' character with a space (' ')
		str[found] = ' ';

		// Find the next occurrence of '\t', starting from the next position
		found = str.find('\t', found + 1);
	}
}

SimpleParsedFile parseBasicFile(const std::string& path, bool verbose, SetSectionFunction setSection, vector<char> ignores = {';', '#'}, char delimiter = ' ')
{
	std::ifstream file;
	file.open(path);
	if (!file.is_open() || file.fail()) {
        throw std::runtime_error(std::format("Failed to open file {}\n", path).c_str());
    }

	SimpleParsedFile parsedfile;

	string current_section = "none";

	int ignore_cnt = 0;

	int skipCnt = 0;

	// Forward declaring for optimization reasons
	string line{}, word{};
	while (getline(file, line)) {

		if (skipCnt > 0) {
			skipCnt--;
			continue;
		}

		replaceTabs(line);

		vector<string> row;
		stringstream ss(line);
		while (getline(ss, word, delimiter)) {
			if (!word.empty()) {
				if (ignoreRow(ignores, word)) {
					break;
				}
				row.push_back(word);
			}
			
		}

		if (row.empty()) { continue; }	// This case happens when a line contains 1 or more spaces, but no words. Space are not regarded as comments, since the separate entries in a line

		bool new_section = setSection(row, current_section, skipCnt);
		if (new_section) { continue; }

		parsedfile.rows.push_back({ current_section, row });
	}

	if (verbose) {
		std::cout << path;
		printf("\n\t%zu rows read. %d rows ignored\n", parsedfile.rows.size(), ignore_cnt);
	}	
	
	return parsedfile;
}
SimpleParsedFile Filehandler::parseItpFile(const std::string& path, bool verbose) {
	assert(path.substr(path.length() - 4) == ".itp");

	SetSectionFunction setSectionFn = [](const std::vector<string>& row, string& current_section, int&) -> bool {
		if (row.size() == 3 && row[0][0] == '[') {

			// Need to handle a special case, because some fuckwits used the same keyword twice - straight to jail!
			if (current_section == "dihedraltypes" && row[1] == "dihedraltypes") {	// Workaround for itp files
				current_section = "improperdihedraltypes";
			}
			else {
				current_section = row[1];
			}

			return true;
		}
		return false;
	};	

	return parseBasicFile(path, verbose, setSectionFn);
}

SimpleParsedFile Filehandler::parseTopFile(const std::string& path, bool verbose)
{
	assert(path.substr(path.length() - 4) == ".top" || path.substr(path.length() - 4) == ".itp");

	SetSectionFunction setSectionFn = [](const std::vector<string>& row, string& current_section, int&) -> bool {
		if (row.size() == 3 && row[0][0] == '[') {

			// Need to handle a special case, because some fuckwits used the same keyword twice - straight to jail!
			if (current_section == "dihedrals" && row[1] == "dihedrals") {	// Workaround for itp files
				current_section = "improperdihedrals";
			}
			else {
				current_section = row[1];
			}

			return true;
		}
		return false;
	};

	std::ifstream file;
	file.open(path);
	if (!file.is_open() || file.fail()) {
		throw std::runtime_error(std::format("Failed to open file {}\n", path).c_str());
	}

	SimpleParsedFile parsedfile;

	string current_section = "none";

	int ignore_cnt = 0;

	int skipCnt = 0;

	// Forward declaring for optimization reasons
	string line{}, word{};
	while (getline(file, line)) {

		if (skipCnt > 0) {
			skipCnt--;
			continue;
		}

		replaceTabs(line);

		vector<char> ignores = { ';', '#' };
		char delimiter = ' ';

		vector<string> row;
		stringstream ss(line);
		while (getline(ss, word, delimiter)) {
			if (!word.empty()) {
				
				if (ignoreRow(ignores, word)) {
					if (word[0] == ';' && current_section == "atoms") {}	// Do nothing, we need these comments to distinct between residues.. i hate these people :(							
					else break;	// Normal case, just skip the comments
				}
				row.push_back(word);
			}

		}

		if (row.empty()) { continue; }	// This case happens when a line contains 1 or more spaces, but no words. Space are not regarded as comments, since the separate entries in a line

		const bool new_section = setSectionFn(row, current_section, skipCnt);
		if (new_section) { continue; }

		parsedfile.rows.push_back({ current_section, row });
	}

	if (verbose) {
		std::cout << path;
		printf("\n\t%zu rows read. %d rows ignored\n", parsedfile.rows.size(), ignore_cnt);
	}

	return parsedfile;
}

SimpleParsedFile Filehandler::parseLffFile(const std::string& path, bool verbose)
{
	assert(path.substr(path.length() - 4) == ".lff");

	SetSectionFunction setSectionFn = [](const std::vector<string>& row, string& current_section, int&) -> bool {
		if (row.size() == 2 && row[0][0] == '#') {
			current_section = row[1];
			return true;
		}
		return false;
	};

	return parseBasicFile(path, verbose, setSectionFn, {'/'}, ' ');
}

SimpleParsedFile Filehandler::parsePrmFile(const std::string& path, bool verbose)
{
	assert(path.substr(path.length() - 4) == ".prm");

	std::vector<std::string> keywords{ "ATOMS", "BONDS", "ANGLES", "DIHEDRALS", "IMPROPER", "NONBONDED", "NBFIX", "CMAP", "END", "HBOND"};	// I despise the inventor of this fkin "system"

	SetSectionFunction setSectionFn = [&](const std::vector<string>& row, string& current_section, int& skipCnt) -> bool {
		if (row[0] == "NONBONDED") {
			skipCnt = 1;	// special case with a line that is not commented - STRAIGHT TO JAIL!
		}

		if (std::find(keywords.begin(), keywords.end(), row[0]) != keywords.end()) {
			current_section = row[0];
			return true;
		}
		return false;
	};

	return parseBasicFile(path, verbose, setSectionFn, { '!' }, ' ');
}

SimpleParsedFile Filehandler::parseGroFile(const std::string& path, bool verbose)
{
	assert(path.substr(path.length() - 4) == ".gro");

	std::ifstream file;
	file.open(path);
	if (!file.is_open() || file.fail()) {
        throw std::runtime_error(std::format("Failed to open file {}\n", path).c_str());
    }


	SimpleParsedFile parsedfile;


	int64_t skipCnt = 2;	// First 2 lines are title and atom count

	const std::vector<int>entry_lengths{ 5, 5, 5, 5, 8, 8, 8 };
	const int min_chars = 5 + 5 + 5 + 5 + 8 + 8 + 8;

	// Forward declaring for optimization reasons
	string line{}, word{};
	while (getline(file, line)) {

		if (skipCnt > 0) {

			if (skipCnt == 2) {
				// 1st line is title
				parsedfile.rows.push_back({ "title" , { line } });
			}
			if (skipCnt == 1) {
				// 2nd line is atom count
				parsedfile.rows.reserve(std::stoi(line));
				parsedfile.rows.push_back({ "n_atoms", {line} });
			}
			
			skipCnt--;
			continue;	// Optim: here by reserving the size at line #2 in the parsedfile.rows
		}

		std::vector<string> row;

		if (line.length() >= min_chars) {
			// Normal case
			//int position = 0;
			//for (auto entry_len : entry_lengths) {
			//	std::string entry = line.substr(position, entry_len);
			//	row.push_back(entry);
			//	position += entry_len;
			//}
			//parsedfile.rows.push_back({ "atoms", row });
			parsedfile.rows.push_back({ "atoms", { line } });

		}
		else if (line.size()> 0){
			// Second last line, with box size
			vector<string> row;
			stringstream ss(line);
			while (std::getline(ss, word, ' ')) {
				if (!word.empty()) {				
					row.push_back(word);
				}
			}
			parsedfile.rows.push_back({ "box_size", row });
		}
	}
	return parsedfile;
}

void Filehandler::createDefaultSimFilesIfNotAvailable(const std::string& dir, float boxsize_nm) {
	const string simparams_path = dir + "/sim_params.txt";
	if (!Filehandler::fileExists(simparams_path)) {	// TODO: Make this string a default-constant somewhere
		const string contents = "";
		std::ofstream file(simparams_path);
		file << contents;
		file.close();
	}

	const string gro_path = dir + "/conf.gro";
	if (!Filehandler::fileExists(gro_path)) {	// TODO: Make this string a default-constant somewhere
		const string boxsize_str = std::to_string(boxsize_nm) + " ";	// add space between dims
		const string contents = " \n0\n\t\t"+ boxsize_str + boxsize_str + boxsize_str;	// Title, n_atoms, box dimensions
		std::ofstream file(gro_path);
		file << contents;
		file.close();
	}

	const string top_path = dir + "/topol.top";
	if (!Filehandler::fileExists(top_path)) {	// TODO: Make this string a default-constant somewhere
		const string contents = "[ moleculetype ]\n[ atoms ]\n[ bonds ]\n[ angles ]\n[ dihedrals ]\n[dihedrals]\n";
		std::ofstream file(top_path);
		file << contents;
		file.close();
	}
}
