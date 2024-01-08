#pragma once

#include "LimaTypes.cuh"
#include "Simulation.cuh"
#include "SimulationDevice.cuh"
#include "EngineUtils.cuh"
#include "BoundaryCondition.cuh"
#include <chrono>
//#include <thread>
//#include <mutex>
//#include <functional>

template <typename BoundaryCondition>
__global__ void updateCompoundNlistsKernel(SimulationDevice* sim_dev);

template <typename BoundaryCondition>
__global__ void updateBlockgridKernel(SimulationDevice* sim_dev);

const int nthreads_in_blockgridkernel = 128;


namespace NeighborLists {


	void updateNlists(SimulationDevice* sim_dev, BoundaryConditionSelect bc, int n_compounds, int& timing);
};
