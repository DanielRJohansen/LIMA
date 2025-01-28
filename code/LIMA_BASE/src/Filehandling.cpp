#include "Filehandling.h"

#include <assert.h>
#include <algorithm>
#include <functional>
#include <array>
#include <fstream>

#include <format>
#include <cctype>
#include <vector>



namespace fs = std::filesystem;
using std::string;


/// <summary>
/// Returns true for a new section. A new section mean the current line (title) is skipped
/// The function can also modify the skipCnt ref so the following x lines are skipped
/// </summary>
using SetSectionFunction = std::function<bool(const std::vector<string>& row, string& section, int& skipCnt)>;


void FileUtils::removeWhitespace(std::string& str) {
	str.erase(std::remove_if(str.begin(), str.end(), 
		[](unsigned char ch) {return std::isspace(static_cast<int>(ch));}),
		str.end());
}

bool FileUtils::firstNonspaceCharIs(const std::string& str, char query) {
	auto first_non_space = std::find_if(str.begin(), str.end(), [](unsigned char ch) {
		return !std::isspace(static_cast<int>(ch));
	});

	return (first_non_space != str.end() && *first_non_space == query);
}

// Reads "key=value" pairs from a file. Disregards all comments (#)
std::unordered_map<std::string, std::string> FileUtils::parseINIFile(const std::string& path, bool forceLowercase) {
	std::ifstream file(path);
	if (!file.is_open()) {
		throw std::runtime_error(std::format("Failed to open file {}\n", path));
	}

	auto ToLowercase = [](std::string& str) {
		std::transform(str.begin(), str.end(), str.begin(), ::tolower);
	};
	

	std::unordered_map<std::string, std::string> dict;
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

			if (forceLowercase) {
				ToLowercase(key);
				ToLowercase(value);
			}

			dict[key] = value;
		}
	}
	return dict;
}

std::string_view FileUtils::ExtractBetweenQuotemarks(const std::string& input) {
	return std::count(input.begin(), input.end(), '"') == 2
		? std::string_view(input.c_str() + input.find('"') + 1, input.find('"', input.find('"') + 1) - input.find('"') - 1)
		: std::string_view{};
}


fs::path FileUtils::GetLimaDir() {
#ifdef __linux__
	return {"/usr/share/LIMA"};
#else
	return { R"(C:\Users\Daniel\git_repo\LIMA)" };
#endif
}

std::vector<std::array<fs::path, 2>> FileUtils::GetAllGroItpFilepairsInDir(const fs::path& dir) {
	std::vector<std::array<fs::path, 2>> pairs;

	for (const auto& entry : fs::directory_iterator(dir)) {
		if (entry.path().extension() == ".gro") {
			std::string base_name = entry.path().stem().string();
			fs::path itp_file = dir / (base_name + ".itp");
			if (fs::exists(itp_file)) {
				pairs.emplace_back(std::array<fs::path, 2>{ entry.path(), itp_file });
			}
		}
	}

	return pairs;
}

void FileUtils::SkipIfdefBlock(std::ifstream& file) {
	std::string line;
	while (std::getline(file, line)) {
		// Trim leading whitespace in-place
		auto pos = line.find_first_not_of(" \t");
		if (pos == std::string::npos) continue; // Skip empty or whitespace-only lines

		// Check for preprocessor directives directly
		if (line.compare(pos, 5, "#else") == 0 || line.compare(pos, 6, "#endif") == 0) {
			return;
		}
	}
	throw std::runtime_error("Failed to find #endif in file\n");
}
bool FileUtils::ChecklineForIfdefAndSkipIfFound(std::ifstream& file, const std::string& line, const std::unordered_set<std::string>& defines) {
	auto pos = line.find_first_not_of(" ");
	if (pos == std::string::npos) return false; // If no non-space character, skip

	if (line.size() >= pos + 6 && line.substr(pos, 6) == "#ifdef") {
		auto keyword = line.substr(pos + 6);
		if (defines.find(keyword) == defines.end()) {
			SkipIfdefBlock(file);
		}
		return true;
	}
	else if (line.size() >= pos + 7 && line.substr(pos, 7) == "#ifndef") {
		auto keyword = line.substr(pos + 7);
		if (defines.find(keyword) != defines.end()) {
			SkipIfdefBlock(file);
		}
		return true;
	}
	return false;
}

std::optional<std::string> FileUtils::ChechlineForDefine(const std::string& line) {
	auto pos = line.find_first_not_of(" ");

	if (pos != std::string::npos) {
		if (line.size() >= pos + 7 && line.substr(pos, 7) == "#define") {
			std::string define = line.substr(pos + 7);
			removeWhitespace(define);
			return define;
		}
	}
	return std::nullopt;
}

std::vector<Float3> FileUtils::ReadCsvAsVectorOfFloat3(const fs::path& path) {
	std::vector<Float3> result;
	std::ifstream file(path);

	if (!file.is_open()) {
		throw std::runtime_error("Could not open file: " + path.string());
	}

	std::string line;
	while (std::getline(file, line)) {
		if (line.empty())
			continue;

		std::istringstream lineStream(line);
		std::string value;
		Float3 temp;

		if (std::getline(lineStream, value, ',')) {
			temp.x = std::stof(value);
		}
		if (std::getline(lineStream, value, ',')) {
			temp.y = std::stof(value);
		}
		if (std::getline(lineStream, value, ',')) {
			temp.z = std::stof(value);
		}

		result.push_back(temp);
	}

	return result;
}

//void FileUtils::SkipIfdefBlock(std::ifstream& file) {
//	std::string line;
//	while (getline(file, line)) {
//		if (line.size() > 5 && line.substr(0, 7) == "#endif") {			
//			return;
//		}
//	}
//
//	throw std::runtime_error(std::format("Failed to find #endif in file\n"));
//}
//
//
//bool FileUtils::ChecklineForIfdefAndSkipIfFound(std::ifstream& file, const std::string& line, const std::unordered_map<std::string>& defines) {
//	if (line.size() > 5 && line.substr(0, 6) == "#ifdef") {
//		SkipIfdefBlock(file);
//		return true;
//	}
//	if (line.size() > 6 && line.substr(0, 7) == "#ifndef") {
//		SkipIfdefBlock(file);
//		return true;
//	}
//	return false;
//}