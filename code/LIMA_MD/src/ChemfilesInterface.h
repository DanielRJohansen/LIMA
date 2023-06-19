#pragma once

#include "LIMA_BASE/include/LimaTypes.cuh"
#include "LIMA_BASE/include/Simulation.cuh"

struct TrrFile {
	static void dumpToFile(const Simulation* sim, const std::string& path);
};


