#pragma once

#include <string>
#include <algorithm>

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