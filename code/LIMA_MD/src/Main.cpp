#include <iostream>

#include "Environment.h"


int main() {
	Environment env{"", EnvMode::Full};
	std::printf("Hello world lima 3");
	return 0;
}
//
//int main() {
//	std::printf("Program starting...\n");
//	LIMA_Print::setScreenSize();		// TODO: fix this for linux
//
//	//Environment::prepFF("conf_test.gro", "topol_test.top");
//
//	const string conf_path = "C:\\PROJECTS\\Quantom\\Simulation\\Molecule\\conf.gro";
//	const string topol_path = "C:\\PROJECTS\\Quantom\\Simulation\\Molecule\\topol.top";
//
//	//Environment::prepFF(conf_path, topol_path);
//
//
//	std::string work_folder = "C:\\PROJECTS\\Quantom\\Simulation\\T4Lysozyme\\";
//
//	Environment env;
//	env.CreateSimulation(conf_path, topol_path, work_folder);
//	env.run();
//
//	//Env.makeVirtualTrajectory("D:\\Quantom\\trajectory.csv", "D:\\Quantom\\waterforce.csv");
//	//Env.renderTrajectory("D:\\Quantom\\virtrj.csv");
//	//Env.renderTrajectory("D:\\Quantom\\trajectory.csv");
//	return 0;
//}

