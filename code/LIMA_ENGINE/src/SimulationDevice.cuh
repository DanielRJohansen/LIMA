#pragma once

#include "Bodies.cuh"
#include "Simulation.cuh"
#include "EngineBodies.cuh"
#include "ChargeOcttree.cuh"

using namespace SCA;

/// <summary>
/// All members of this will only ever exist on device. Immediately after creating this class
/// it must also be moved to device.
/// </summary>
struct SimulationDevice {
	SimulationDevice(const SimulationDevice&) = delete;

	SimulationDevice(const SimParams& params_host, std::unique_ptr<Box> box_host)
	{
		// Allocate structures for keeping track of solvents and compounds
		//cudaMallocManaged(&compound_grid, sizeof(CompoundGrid));		
		compound_grid = CompoundGrid::MallocOnDevice(box_host->boxparams.boxSize);
		cudaMallocManaged(&compound_neighborlists, sizeof(NeighborList) * MAX_COMPOUNDS);
		cudaMallocManaged(&transfermodule_array, sizeof(SolventBlockTransfermodule) * SolventBlocksCircularQueue::blocks_per_grid);

		{
			SimSignals temp{};
			genericCopyToDevice(temp, &signals, 1);
		}
		
		genericCopyToDevice(params_host, &params, 1);

		databuffers = new DatabuffersDevice(box_host->boxparams.total_particles_upperbound, box_host->boxparams.n_compounds, params_host.data_logging_interval);
		databuffers = genericMoveToDevice(databuffers, 1);

		box_host->moveToDevice();
		cudaMallocManaged(&box, sizeof(Box));
		cudaMemcpy(box, box_host.get(), sizeof(Box), cudaMemcpyHostToDevice);

		ChargeOctTree charge_octtree_host;
		cudaMallocManaged(&charge_octtree, sizeof(ChargeOctTree));
		cudaMemcpy(charge_octtree, &charge_octtree_host, sizeof(ChargeOctTree), cudaMemcpyHostToDevice);

		box_host->owns_members = false;
		box_host->is_on_device = false; // because moveToDevice sets it to true before transferring.
		box_host.reset();
	}

	// Recursively free members. Use cudaFree on *this immediately after
	void deleteMembers() {
		box->deleteMembers();
		cudaFree(box);

		databuffers->freeMembers();
		cudaFree(databuffers);

		cudaFree(params);

		CompoundGrid::Free(compound_grid);
		cudaFree(compound_neighborlists);

		cudaFree(charge_octtree);
	}

	// Compounds signal where they are on a grid, handled by NLists. Used by solvents to load nearby compounds.
	CompoundGrid* compound_grid = nullptr;

	// Compounds can see which compounds are near them
	NeighborList* compound_neighborlists = nullptr;

	// Module used to move solvents to a new block, in parallel
	SolventBlockTransfermodule* transfermodule_array = nullptr;

	SimParams* params;
	SimSignals* signals;
	Box* box;
	DatabuffersDevice* databuffers;

	ChargeOctTree* charge_octtree;
};