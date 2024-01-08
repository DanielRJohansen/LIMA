#pragma once

#include "Bodies.cuh"
#include "BoundaryCondition.cuh"	// TODO: Remove
#include "Simulation.cuh"


struct SimulationDevice {
	SimulationDevice(const SimulationDevice&) = delete;

	SimulationDevice(const SimParams& params_host, std::unique_ptr<Box> box_host) {
		genericCopyToDevice(params_host, &params, 1);

		databuffers = new DatabuffersDevice(box_host->boxparams.total_particles_upperbound, box_host->boxparams.n_compounds);
		databuffers = genericMoveToDevice(databuffers, 1);

		box_host->moveToDevice();
		cudaMallocManaged(&box, sizeof(Box));
		cudaMemcpy(box, box_host.get(), sizeof(Box), cudaMemcpyHostToDevice);

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
	}
	SimParams* params;
	Box* box;
	DatabuffersDevice* databuffers;
};