#pragma once

#include <string>
#include <algorithm>

namespace CmdLineUtils {
	char* getCmdOption(char** begin, char** end, const std::string& option)
	{
		char** itr = std::find(begin, end, option);
		if (itr != end && std::next(itr) != end) {
			return *std::next(itr);
		}
		return nullptr;
	}

	bool cmdOptionExists(char** begin, char** end, const std::string& option)
	{
		return std::find(begin, end, option) != end;
	}

	std::string ToLowercase(char* str)
	{
		std::string s(str);
		std::transform(s.begin(), s.end(), s.begin(), ::tolower);
		return s;
	}
}