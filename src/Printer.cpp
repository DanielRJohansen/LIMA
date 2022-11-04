#include "Printer.h"

#include <math.h>
#include <algorithm>


//using namespace LIMA_Printer;

void addSpaces(std::string& str, int n_spaces) {
	for (int i = 0; i < n_spaces; i++)
		str += " ";
}

std::string formatValue(int value)  {
	return std::to_string(value);
}

std::string formatValue(float value) {
	std::string val_s = std::to_string(value);

	int decimals = 6 - log10(value);	
	std::string val_s_rounded = val_s.substr(0, val_s.find(".") + decimals);
	
	return val_s_rounded;
}

std::string formatValue(std::variant<int, float>& value) {
	if (holds_alternative<int>(value)) { return formatValue(std::get<int>(value)); }
	if (holds_alternative<int>(value)) { return formatValue(std::get<float>(value)); }
	printf("Bad value type\n");
	exit(0);
}

void addRightadjustedStringToString(std::string& main_string, std::string& str, const int chars_per_elem) {
	addSpaces(main_string, chars_per_elem - str.size());
	main_string += str;
}

//
//
//void LIMA_Printer::printNameValuePairs(std::vector < std::pair < std::string, std::variant<int, float>>> matrix) {
//	const int chars_per_elem = default_width / matrix.size();
//
//	std::vector<std::string> lines{"", ""};
//
//	for (auto pair : matrix) {
//
//		addRightadjustedStringToString(lines[0], pair.first, chars_per_elem);
//
//		std::string value_as_string = formatValue(pair.second);
//		addRightadjustedStringToString(lines[1], value_as_string, chars_per_elem);
//	}
//}

//namespace LIMA_Printer {
//	void test(std::variant<int, float> val) {}
//}
