#pragma once
#include <vector>
#include <iostream>
#include <variant>
#include <string>




class LIMA_Printer {
public:

	//void test(std::variant<int, float> val) {}

	//void printNameValuePairs(std::vector < std::pair < std::string, std::variant<int, float>>> matrix);


	void printH1(std::string);
	void printH2(std::string);






	// sizes in chars
	const int default_height = 25;
	const int default_width = 80;


};