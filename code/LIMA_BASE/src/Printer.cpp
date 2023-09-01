#include "Printer.h"

#include <math.h>

#ifndef __linux__
#include <Windows.h>
#endif

//using namespace LIMA_Printer;


void addMultipleChars(std::string& str, int n_spaces, char c = ' ') {
	for (int i = 0; i < n_spaces; i++)
		str += string{ c };
}


//std::string LIMA_Printer::formatValue(int value)  {
//	return std::to_string(value);
//}

std::string LIMA_Printer::formatValue(double value) {
	return formatValue(static_cast<float>(value));
}

std::string LIMA_Printer::formatValue(float value)  {
	std::string val_s = std::to_string(value);

	int decimals = 6 - static_cast<int>(log10(value));
	std::string val_s_rounded = val_s.substr(0, val_s.find(".") + decimals);

	return val_s_rounded;
}



void LIMA_Printer::addRightadjustedStringToString(std::string& main_string, const std::string& str) {
	addMultipleChars(main_string, static_cast<int>(chars_per_elem - str.size()));
	main_string += str;
}

void LIMA_Printer::printTableRow(const std::vector<string>& row) {
	string str = "";
	for (const auto& elem : row) {
		addRightadjustedStringToString(str, elem);
	}
	std::cout << str << "\n\n";
}

void LIMA_Printer::printTableRow(string s, std::vector<float> data) {
	std::vector<string> row;
	row.push_back(s);
	for (auto elem : data) {
		row.push_back(formatValue(elem));
	}
	printTableRow(row);
}






// Namespace
	// sizes in chars
static const int default_height = 60;
static const int default_width = 120;
static const int chars_per_elem = default_width / 6;

void LIMA_Print::setScreenSize()
{
#ifndef __linux__
	//HWND hwnd = GetConsoleWindow();
	//if (hwnd != NULL) { MoveWindow(hwnd, 0, 0, default_width, default_height, TRUE); }
#endif
}

void LIMA_Print::printH(std::string str, char c, bool ls, bool ts) {
	string out_str = "";
	if (ls) out_str += "\n";
	if (str == "") { 
		addMultipleChars(out_str, default_width, c); 
	}
	else {
		const int n_fillchars_total = static_cast<int>(default_width - str.size() - 2);
		const int extra_char_left = n_fillchars_total % 2;

		addMultipleChars(out_str, n_fillchars_total / 2 + extra_char_left, c);
		out_str += " ";
		out_str += str;
		out_str += " ";
		addMultipleChars(out_str, n_fillchars_total / 2, c);
	}
	out_str += "\n";
	if (ts) out_str += "\n";
	std::cout << out_str;
}

void LIMA_Print::printH1(std::string str, bool ls, bool ts) {
	printH(str, '#', ls, ts);
}
void LIMA_Print::printH2(std::string str, bool ls, bool ts) {
	printH(str, '-', ls, ts);
}
