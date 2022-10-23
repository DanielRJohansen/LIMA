#include <iostream>
#include <math.h>

#include "Environment.h"
#include "Engine.cuh"

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




int main() {

	int a = 0b1011;
	int b = 0b1 << 0;
	int c = a & b;
	printf("%d", c);
	exit(0);
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

