#pragma once
#include <vector>
#include <iostream>

#include <string>



using string = std::string;

namespace LIMA_Print {
	void setScreenSize();
	void printH(std::string, char c, bool leading_space, bool trailing_space);
	void printH1(std::string = "", bool ls = false, bool ts = true);
	void printH2(std::string = "", bool ls = false, bool ts = true);

	template<typename T>
	void printMatlabVec(std::string name, const std::vector<T>& vec) {
		std::cout << name << " = [";
		for (const auto& elem : vec) { 
			std::cout << elem << " "; 
		}
		std::cout << "];\n";
	}
}

class LIMA_Printer {
public:

	//void test(std::variant<int, float> val) {}

	//void printNameValuePairs(std::vector < std::pair < std::string, std::variant<int, float>>> matrix);

	template <typename T>
	static void doThing(string str, T val, std::vector<string>& buffer) {
//		std::vector<std::string> lines{ "", "" };
		addRightadjustedStringToString(buffer[0], str);
		std::string value_as_string = formatValue(val);

		addRightadjustedStringToString(buffer[1], value_as_string);

	}

	//template void doThing<int>(std::pair<string, T> p1, std::vector<string>& buffer);

	template <typename T>
	static void printNameValuePairs(string s1, T v1, std::vector<string> buffer = {"", ""}) {
		doThing(s1, v1, buffer);
		std::cout << buffer[0] << "\n" << buffer[1] << "\n";
	}

	template <typename T1, typename T2>
	static void printNameValuePairs(string s1, T1 v1, string s2, T2 v2, std::vector<string> buffer = { "", "" }) {
		doThing(s1, v1, buffer);
		printNameValuePairs<T2>(s2, v2, buffer);
	}

	template <typename T1, typename T2, typename T3>
	static void printNameValuePairs(string s1, T1 v1, string s2, T2 v2, string s3, T3 v3, std::vector<string> buffer = { "", "" }) {
		doThing(s1, v1, buffer);
		printNameValuePairs<T2, T3>(s2, v2, s3, v3, buffer);
	}

	template <typename T1, typename T2, typename T3, typename T4>
	static void printNameValuePairs(string s1, T1 v1, string s2, T2 v2, string s3, T3 v3, string s4, T4 v4, std::vector<string> buffer = { "", "" }) {
		doThing(s1, v1, buffer);
		printNameValuePairs<T2, T3, T4>(s2, v2, s3, v3, s4, v4, buffer);
	}

	template <typename T1, typename T2, typename T3, typename T4, typename T5>
	static void printNameValuePairs(string s1, T1 v1, string s2, T2 v2, string s3, T3 v3, string s4, T4 v4, string s5, T5 v5, std::vector<string> buffer = { "", "" }) {
		doThing(s1, v1, buffer);
		printNameValuePairs<T2, T3, T4, T5>(s2, v2, s3, v3, s4, v4, s5, v5, buffer);
	}

	static void printTableRow(const std::vector<string>& row);
	static void printTableRow(string s, std::vector<float> data);


private:
	
	// I have to do this because CUDA is being a bitch..
	//static std::string formatValue(int value);
	template <typename T>
	static std::string formatValue(T value) {
		return std::to_string(value);
	}
	static std::string formatValue(float value);
	static std::string formatValue(double value);


	static void addRightadjustedStringToString(std::string& main_string, const std::string& str);


	// sizes in chars
	static const int default_height = 60;
	static const int default_width = 120;
	static const int chars_per_elem = default_width / 6;
};

