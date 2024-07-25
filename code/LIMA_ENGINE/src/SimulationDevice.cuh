#pragma once

// TOdo move impl to .cpp file

#include "Bodies.cuh"
#include "Simulation.cuh"
#include "EngineBodies.cuh"
#include "ChargeOcttree.cuh"

//using namespace SCA;

/// <summary>
/// All members of this will only ever exist on device. Immediately after creating this class
/// it must also be moved to device.
/// </summary>
struct SimulationDevice {
	SimulationDevice(const SimulationDevice&) = delete;

	SimulationDevice(const SimParams& params_host, std::unique_ptr<Box> box_host)
	{
		// Allocate structures for keeping track of solvents and compounds
		compound_grid = BoxGrid::MallocOnDevice<CompoundGridNode>(box_host->boxparams.boxSize);
		cudaMallocManaged(&compound_neighborlists, sizeof(NeighborList) * MAX_COMPOUNDS);
		
		cudaMallocManaged(&transfermodule_array, sizeof(SolventBlockTransfermodule) * BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(box_host->boxparams.boxSize)));

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

		//ChargeOctTree charge_octtree_host;
		//cudaMallocManaged(&charge_octtree, sizeof(ChargeOctTree));
		//cudaMemcpy(charge_octtree, &charge_octtree_host, sizeof(ChargeOctTree), cudaMemcpyHostToDevice);


		chargeGrid = BoxGrid::MallocOnDevice<Electrostatics::ChargeNode>(box_host->boxparams.boxSize);
		chargeGridChargeSums = BoxGrid::MallocOnDevice<float>(box_host->boxparams.boxSize);
		chargeGridOutputForceAndPot = BoxGrid::MallocOnDevice<ForceAndPotential>(box_host->boxparams.boxSize);

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

		cudaFree(compound_grid);
		cudaFree(compound_neighborlists);

		//cudaFree(charge_octtree);
		cudaFree(chargeGrid);
		cudaFree(chargeGridChargeSums);
		cudaFree(chargeGridOutputForceAndPot);
	}

	// Compounds signal where they are on a grid, handled by NLists. Used by solvents to load nearby compounds.
	CompoundGridNode* compound_grid = nullptr;

	// Compounds can see which compounds are near them
	NeighborList* compound_neighborlists = nullptr;

	// Module used to move solvents to a new block, in parallel
	SolventBlockTransfermodule* transfermodule_array = nullptr;

	SimParams* params;
	SimSignals* signals;
	Box* box;
	DatabuffersDevice* databuffers;

	//ChargeOctTree* charge_octtree;
	Electrostatics::ChargeNode* chargeGrid = nullptr;
	float* chargeGridChargeSums = nullptr;

	// potE should be divided equally between all the particles in the node
	ForceAndPotential* chargeGridOutputForceAndPot = nullptr; // {Float3 force [1/l N/mol], float potE [J/mol]}

};