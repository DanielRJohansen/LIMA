#include "MDCalibrations.h"
#include "Environment.h"

#include "LimaTypes.cuh"

#include <iostream>
#include <string>




bool doPoolCalibration(Environment& env) {
	const std::string work_folder = "C:/PROJECTS/Quantom/Simulation/Pool/";
	const std::string conf = work_folder + "molecule/conf.gro";
	const std::string topol = work_folder + "molecule/topol.top";

	//env.CreateSimulation(conf, topol, work_folder);
}




int main() {
	Environment env;

	//doPoolCalibration(env);

	Float3 a{ 0.f };


	return 0;
}