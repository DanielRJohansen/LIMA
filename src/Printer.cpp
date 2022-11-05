#include "Printer.h"

#include <math.h>



//using namespace LIMA_Printer;

void addSpaces(std::string& str, int n_spaces) {
	for (int i = 0; i < n_spaces; i++)
		str += " ";
}

std::string LIMA_Printer::formatValue(int value)  {
	return std::to_string(value);
}

std::string LIMA_Printer::formatValue(double value) {
	formatValue(static_cast<float>(value));
}

std::string LIMA_Printer::formatValue(float value)  {
	std::string val_s = std::to_string(value);

	int decimals = 6 - log10(value);
	std::string val_s_rounded = val_s.substr(0, val_s.find(".") + decimals);

	return val_s_rounded;
}



void LIMA_Printer::addRightadjustedStringToString(std::string& main_string, std::string& str) {
	addSpaces(main_string, chars_per_elem - str.size());
	main_string += str;
}
