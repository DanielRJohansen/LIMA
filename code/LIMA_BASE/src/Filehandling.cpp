#include "Filehandling.h"

#include <assert.h>
#include <algorithm>
#include <functional>
#include <map>
#include <array>
#include <fstream>

#include <format>
#include <cctype>
#include <vector>



#ifdef __linux__
#include <unistd.h>
#else
#include <windows.h>
#endif






namespace fs = std::filesystem;
using std::string;


/// <summary>
/// Returns true for a new section. A new section mean the current line (title) is skipped
/// The function can also modify the skipCnt ref so the following x lines are skipped
/// </summary>
using SetSectionFunction = std::function<bool(const std::vector<string>& row, string& section, int& skipCnt)>;

bool ignoreRow(const std::vector<char>& ignores, const string& line) {
	if (line.length() == 0)
		return true;
	for (auto& c : ignores) {
		if (line[0] == c)
			return true;
	}
	return false;
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

// Reads "key=value" pairs from a file. Disregards all comments (#)
std::map<std::string, std::string> Filehandler::parseINIFile(const std::string& path) {
	std::ifstream file(path);
	if (!file.is_open()) {
		throw std::runtime_error(std::format("Failed to open file {}\n", path));
	}

	std::map<std::string, std::string> dict;
	std::string line;
	while (getline(file, line)) {
		// Sanitize line by removing everything after a '#' character
		size_t commentPos = line.find('#');
		if (commentPos != std::string::npos) {
			line = line.substr(0, commentPos);
		}

		std::stringstream ss(line);
		std::string key, value;
		if (getline(ss, key, '=') && getline(ss, value)) {
			// Remove spaces from both key and value
			key.erase(remove_if(key.begin(), key.end(), ::isspace), key.end());
			value.erase(remove_if(value.begin(), value.end(), ::isspace), value.end());
			dict[key] = value;
		}
	}
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

SimpleParsedFile parseBasicFile(const std::string& path, bool verbose, SetSectionFunction setSection, std::vector<char> ignores = {';', '#'}, char delimiter = ' ')
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

		std::vector<string> row;
		std::stringstream ss(line);
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


void Filehandler::createDefaultSimFilesIfNotAvailable(const std::string& dir, float boxsize_nm) {
	const string simparams_path = dir + "/sim_params.txt";
	if (!std::filesystem::exists(simparams_path)) {	// TODO: Make this string a default-constant somewhere
		const string contents = "";
		std::ofstream file(simparams_path);
		file << contents;
		file.close();
	}

	const string gro_path = dir + "/conf.gro";
	if (!std::filesystem::exists(gro_path)) {	// TODO: Make this string a default-constant somewhere
		const string boxsize_str = std::to_string(boxsize_nm) + " ";	// add space between dims
		const string contents = " \n0\n\t\t"+ boxsize_str + boxsize_str + boxsize_str;	// Title, nAtoms, box dimensions
		std::ofstream file(gro_path);
		file << contents;
		file.close();
	}

	const string top_path = dir + "/topol.top";
	if (!std::filesystem::exists(top_path)) {	// TODO: Make this string a default-constant somewhere
		const string contents = "[ moleculetype ]\n[ atoms ]\n[ bonds ]\n[ angles ]\n[ dihedrals ]\n[dihedrals]\n";
		std::ofstream file(top_path);
		file << contents;
		file.close();
	}
}

fs::path Filehandler::GetLimaDir() {
#ifdef __linux__
	return {"/usr/share/LIMA"};
#else
	const int maxPathLen = 2048;
	char executablePath[maxPathLen];
	GetModuleFileName(NULL, executablePath, maxPathLen);


	// Start from the directory of the executable
	fs::path execPath(executablePath);
	fs::path currentPath = execPath.parent_path();

	// Traverse up the directory tree until "LIMA" directory is found or root is reached
	int count = 0;
	while (currentPath.has_parent_path()) {
		if (currentPath.filename() == "LIMA") {
			return currentPath; // Return the path when "LIMA" is found
		}
		currentPath = currentPath.parent_path(); // Move one level up in the directory tree
		if (count++ > 20)
			break;
	}

	throw std::runtime_error("Could not find LIMA directory");
#endif
}
