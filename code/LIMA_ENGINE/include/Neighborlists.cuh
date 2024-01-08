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

	// Must come before any other kernel()
	void _updateCompoundNlistsKernel(int n_compounds, SimulationDevice* sim_dev);
	void _updateBlockgridKernel(int n_compounds, SimulationDevice* sim_dev);

	template <typename BoundaryCondition>
	static void updateNlists(SimulationDevice* sim_dev, int n_compounds, int& timing)
	{
		const auto t0 = std::chrono::high_resolution_clock::now();

		BoundaryCondition bc;

		{
			_updateCompoundNlistsKernel(n_compounds, sim_dev);
			//const int n_blocks = n_compounds / threads_in_compoundnlist_kernel + 1;
			//updateCompoundNlistsKernel << < n_blocks, threads_in_compoundnlist_kernel >> > (sim_dev);	// Must come before any other kernel()
		}

		cudaDeviceSynchronize();	// The above kernel overwrites the nlists, while the below fills ut the nlists present, so the above must be completed before progressing

#ifdef ENABLE_SOLVENTS
		{
			_updateBlockgridKernel(n_compounds, sim_dev);
			//const int n_blocks = CompoundGrid::blocks_total / nthreads_in_blockgridkernel + 1;
			//updateBlockgridKernel << <n_blocks, nthreads_in_blockgridkernel >> > (sim_dev);
		}
#endif

		cudaDeviceSynchronize();
		const auto t1 = std::chrono::high_resolution_clock::now();
		timing += static_cast<int>(std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count());
	}

	//template <>	void updateNlists<PeriodicBoundaryCondition>
	//(SimulationDevice<PeriodicBoundaryCondition>* sim_dev, int n_compounds, int& timing);


	//template <>	void updateNlists<NoBoundaryCondition>(SimulationDevice<NoBoundaryCondition>* sim_dev, int n_compounds, int& timing);

};
