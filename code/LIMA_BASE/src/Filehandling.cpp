#include "Filehandling.h"

#include <assert.h>
#include <algorithm>
#include <functional>
#include <array>
#include <fstream>
#include <mutex>
#include <optional>

#include <format>
#include <cctype>
#include <vector>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
#include <unistd.h>
#endif


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

bool FileUtils::firstNonspaceCharIs(const std::string_view& str, char query) {
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


namespace {
	fs::path ExecutableDir() {
#ifdef _WIN32
		wchar_t buffer[MAX_PATH];
		const DWORD length = GetModuleFileNameW(nullptr, buffer, MAX_PATH);
		return length > 0 && length < MAX_PATH ? fs::path(buffer).parent_path() : fs::path{};
#else
		std::error_code error;
		const fs::path executable = fs::read_symlink("/proc/self/exe", error);
		return error ? fs::path{} : executable.parent_path();
#endif
	}

	std::optional<fs::path> FindRepositoryRoot(fs::path path) {
		while (!path.empty() && path != path.parent_path()) {
			if (fs::exists(path / "lima_main_dir.txt"))
				return path;
			path = path.parent_path();
		}
		return std::nullopt;
	}

	fs::path FindLimaDir() {
		const fs::path executableDir = ExecutableDir();
		if (!executableDir.empty()) {
			// Windows release: resources next to lima.exe
			if (fs::exists(executableDir / "resources"))
				return executableDir;
			// Linux packages and tarball: bin/lima with share/LIMA/resources
			if (fs::exists(executableDir.parent_path() / "share" / "LIMA" / "resources"))
				return executableDir.parent_path() / "share" / "LIMA";
			// Development builds live inside the repository
			if (auto root = FindRepositoryRoot(executableDir))
				return *root;
		}
		if (auto root = FindRepositoryRoot(fs::current_path()))
			return *root;
#ifdef __linux__
		if (fs::exists("/usr/share/LIMA/resources"))
			return "/usr/share/LIMA";
#endif
		throw std::runtime_error(std::format(
			"Could not find LIMA's resources directory. Searched next to the executable ({}), in ../share/LIMA, and in the parent directories",
			executableDir.string()));
	}
}

fs::path FileUtils::GetLimaDir() {
	// Thread safe, as it is called concurrently by Environment's worker threads
	static const fs::path limaDir = FindLimaDir();
	return limaDir;
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

std::string FileUtils::ReadFileToString(const fs::path& path) {
	std::ifstream file(path, std::ios::binary | std::ios::ate);
	if (!file)
		throw std::runtime_error(std::format("Failed to open file {}\n", path.string()));

	const std::streamsize size = file.tellg();
	std::string contents(static_cast<size_t>(size), '\0');

	file.seekg(0);
	file.read(contents.data(), size);

	if (!file)
		throw std::runtime_error(std::format("Failed to read file {}\n", path.string()));

	return contents;
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