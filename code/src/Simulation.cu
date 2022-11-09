#include "Simulation.cuh"


void Box::moveToDevice() {
	int bytes_total = sizeof(Compound) * n_compounds
		+ sizeof(Solvent) * MAX_SOLVENTS * 2
		+ sizeof(CompoundState) * MAX_COMPOUNDS * 2
		+ sizeof(NeighborList) * (MAX_SOLVENTS + MAX_COMPOUNDS);
	printf("BOX: moving %.2f MB to device\n", (float)bytes_total * 1e-6);

	compounds = genericMoveToDevice(compounds, n_compounds);
	solvents = genericMoveToDevice(solvents, MAX_SOLVENTS);
	bridge_bundle = genericMoveToDevice(bridge_bundle, 1);
	compound_state_array = genericMoveToDevice(compound_state_array, MAX_COMPOUNDS);


	compound_neighborlists = genericMoveToDevice(compound_neighborlists, MAX_COMPOUNDS);
	solvent_neighborlists = genericMoveToDevice(solvent_neighborlists, MAX_SOLVENTS);

	bonded_particles_lut_manager = genericMoveToDevice(bonded_particles_lut_manager, 1);

	cudaDeviceSynchronize();
	printf("Box transferred to device\n");
}


void Simulation::moveToDevice() {
	box = genericMoveToDevice(box, 1);
	cudaDeviceSynchronize();
	if (cudaGetLastError() != cudaSuccess) {
		fprintf(stderr, "Error during Simulation Host->Device transfer\n");
		exit(1);
	}
	printf("Simulation ready for device\n\n");
}

void Simulation::copyBoxVariables() {
	n_compounds = box->n_compounds;
	n_bridges = box->bridge_bundle->n_bridges;


	n_solvents = box->n_solvents;
	blocks_per_solventkernel = (int)ceil((float)n_solvents / (float)THREADS_PER_SOLVENTBLOCK);


	compounds_host = new Compound[n_compounds];
	for (int i = 0; i < n_compounds; i++)
		compounds_host[i] = box->compounds[i];
}