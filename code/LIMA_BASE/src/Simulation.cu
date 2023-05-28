#include "Simulation.cuh"


void Box::moveToDevice() {
	int bytes_total = sizeof(Compound) * n_compounds
		+ sizeof(CompoundState) * MAX_COMPOUNDS * 3
		+ sizeof(NeighborList) * (MAX_SOLVENTS + MAX_COMPOUNDS);
	printf("BOX: moving %.2f MB to device\n", (float)bytes_total * 1e-6);

	compounds = genericMoveToDevice(compounds, n_compounds);
	bridge_bundle = genericMoveToDevice(bridge_bundle, 1);

	coordarray_circular_queue = genericMoveToDevice(coordarray_circular_queue, MAX_COMPOUNDS * STEPS_PER_LOGTRANSFER);
	solventblockgrid_circular_queue = genericMoveToDevice(solventblockgrid_circular_queue, STEPS_PER_SOLVENTBLOCKTRANSFER);


	compound_neighborlists = genericMoveToDevice(compound_neighborlists, MAX_COMPOUNDS);
	solvent_neighborlists = genericMoveToDevice(solvent_neighborlists, MAX_SOLVENTS);

	bonded_particles_lut_manager = genericMoveToDevice(bonded_particles_lut_manager, 1);

	forcefield = genericMoveToDevice(forcefield, 1);

	cudaDeviceSynchronize();
	printf("Box transferred to device\n");
}

Box SimUtils::copyToHost(const Box* box_dev) {
	Box box{};
	cudaMemcpy(&box, box_dev, sizeof(Box), cudaMemcpyDeviceToHost);

	box.compounds = genericCopyToHost(box.compounds, box.n_compounds);


}

Simulation::Simulation(InputSimParams& ip) :
	simparams_host{ SimParamsConst{ip.n_steps, ip.dt} }
{
	box = new Box();
}

Simulation::Simulation(const Box& inputbox, const uint32_t inputbox_current_step, const InputSimParams& ip) :
	simparams_host{ SimParamsConst{ ip.n_steps, ip.dt } }
{
	//box = new Box(inputbox, inputbox_current_step);
}

Simulation::~Simulation() {
	deleteBoxMembers();	// TODO: move to box destructor

}

void Simulation::moveToDevice() {
	box = genericMoveToDevice(box, 1);

	genericCopyToDevice(simparams_host, &simparams_device, 1);

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
		cudaFree(box);
	}
	else {
		delete[] box->compounds;	// TODO: Finish this
		//delete[] box->solvents;
		delete[] box->bridge_bundle;
		delete[] box->compound_neighborlists;
		delete[] box->solvent_neighborlists;
		delete[] box->bonded_particles_lut_manager;
		delete[] box->compounds;
		delete box;
	}
}

void InputSimParams::overloadParams(std::map<std::string, double>& dict) {
	overloadParam(dict, &dt, "dt", FEMTO_TO_LIMA);	// convert [fs] to [ls]
	overloadParam(dict, &n_steps, "n_steps");
}