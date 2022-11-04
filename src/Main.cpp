#include <iostream>
#include <math.h>

#include "Environment.h"
#include "Engine.cuh"

#include <vector>

struct Test {

	int arr[4];
	//const static int size = 14;
	volatile int a = 10;
	bool redu(int b) {
		if (b % 2) {
			printf("Bool %d\n", b % 2);
			a--;
			return 1;
		}
		return 0;
	}
	int get() { return a; }
};


#include "Printer.h"

int main() {

	int a = 0b1011;
	int b = 0b1 << 0;
	int c = a & b;

	std::vector<int> a = { 1, 23, 3 };
	std::vector<std::pair<std::string, int>> b{ {"hey", 2} };
	//LIMA_Printer::printNameValuePairs(b);
	//LIMA_Printer::printNameValuePairs({ {"hey", 2} });


	return 0;

	//Environment::prepFF("conf_test.gro", "topol_test.top");

	string conf_filename = "conf.gro";
	string topol_filename = "topol.top";
	//Environment::prepFF("conf.gro", "topol.top");
	Environment::prepFF(conf_filename, topol_filename);
	std::printf("Program starting...\n");



	Environment Env(conf_filename, topol_filename);
	Env.run();

	//Env.makeVirtualTrajectory("D:\\Quantom\\trajectory.csv", "D:\\Quantom\\waterforce.csv");
	//Env.renderTrajectory("D:\\Quantom\\virtrj.csv");
	//Env.renderTrajectory("D:\\Quantom\\trajectory.csv");
	return 0;
}

