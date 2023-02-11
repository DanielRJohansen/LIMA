#include "Simulation.cuh"


void Box::moveToDevice() {
	int bytes_total = sizeof(Compound) * n_compounds
		//+ sizeof(Solvent) * MAX_SOLVENTS * 2
		+ sizeof(CompoundState) * MAX_COMPOUNDS * 3
		+ sizeof(NeighborList) * (MAX_SOLVENTS + MAX_COMPOUNDS);
	printf("BOX: moving %.2f MB to device\n", (float)bytes_total * 1e-6);

	compounds = genericMoveToDevice(compounds, n_compounds);
	//solvents = genericMoveToDevice(solvents, MAX_SOLVENTS);
	bridge_bundle = genericMoveToDevice(bridge_bundle, 1);

	compound_neighborlists = genericMoveToDevice(compound_neighborlists, MAX_COMPOUNDS);
	solvent_neighborlists = genericMoveToDevice(solvent_neighborlists, MAX_SOLVENTS);

	bonded_particles_lut_manager = genericMoveToDevice(bonded_particles_lut_manager, 1);

	forcefield_device_box = genericMoveToDevice(forcefield_device_box, 1);

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
	box_is_on_device = true;
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

	// Need this variable both on host and device
	//total_particles_upperbound = box->n_compounds * MAX_COMPOUND_PARTICLES + box->n_solvents;											// BAD AMBIGUOUS AND WRONG CONSTANTS	
	total_particles_upperbound = box->n_compounds * MAX_COMPOUND_PARTICLES + SolventBlockGrid::blocks_total * MAX_SOLVENTS_IN_BLOCK;
	box->total_particles_upperbound = total_particles_upperbound;
}

Simulation::Simulation(SimulationParams& sim_params) :
	dt(sim_params.dt),	// convert [ns] to [fs]
	n_steps(sim_params.n_steps)
{
	box = new Box();
}
Simulation::~Simulation() {
	deleteBoxMembers();
	delete box;
}

void Simulation::deleteBoxMembers() {
	if (box_is_on_device) {
		cudaFree(box->compounds);

		cudaFree(box->compound_neighborlists);
		cudaFree(box->solvent_neighborlists);

		//cudaFree(box->solvents);
		//cudaFree(box->solvents_next);

		cudaFree(box->bridge_bundle);		
		cudaFree(box->bonded_particles_lut_manager);

		cudaFree(box->potE_buffer);
		cudaFree(box->traj_buffer);
		cudaFree(box->outdata);
		cudaFree(box->data_GAN);
	}
	else {
		delete[] box->compounds;	// TODO: Finish this
		//delete[] box->solvents;
		delete[] box->bridge_bundle;
		delete[] box->compound_neighborlists;
		delete[] box->solvent_neighborlists;
		delete[] box->bonded_particles_lut_manager;
		delete[] box->compounds;
	}
}

void SimulationParams::overloadParams(std::map<std::string, double>& dict) {
	overloadParam(dict, &dt, "dt", FEMTO_TO_LIMA);	// convert [fs] to [ls]
	overloadParam(dict, &n_steps, "n_steps");
}