#pragma once

#include "LimaTypes.cuh"
#include "Simulation.cuh"

struct TrrFile {
	static void dumpToFile(const Simulation* sim, const std::string& path);
};


