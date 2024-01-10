#pragma once

#include "Simulation.cuh"

class SimulationDevice;

template <typename BoundaryCondition>
__global__ void updateCompoundNlistsKernel(SimulationDevice* sim_dev);

template <typename BoundaryCondition>
__global__ void updateBlockgridKernel(SimulationDevice* sim_dev);



namespace NeighborLists {
	void updateNlists(SimulationDevice*, BoundaryConditionSelect, const BoxParams&, int& timing);
};
