#include "Simulation.cuh"


void Box::moveToDevice() {
	int bytes_total = sizeof(Compound) * n_compounds
		+ sizeof(CompoundState) * MAX_COMPOUNDS * 3
		+ sizeof(NeighborList) * (MAX_SOLVENTS + MAX_COMPOUNDS);
	printf("BOX: moving %.2f MB to device\n", (float)bytes_total * 1e-6);

	compounds = genericMoveToDevice(compounds, MAX_COMPOUNDS);
	bridge_bundle = genericMoveToDevice(bridge_bundle, 1);

	coordarray_circular_queue = genericMoveToDevice(coordarray_circular_queue, MAX_COMPOUNDS * STEPS_PER_LOGTRANSFER);
	solventblockgrid_circular_queue = genericMoveToDevice(solventblockgrid_circular_queue, STEPS_PER_SOLVENTBLOCKTRANSFER);


	compound_neighborlists = genericMoveToDevice(compound_neighborlists, MAX_COMPOUNDS);
	//solvent_neighborlists = genericMoveToDevice(solvent_neighborlists, MAX_SOLVENTS);		// TODO: are we still using these??!

	bonded_particles_lut_manager = genericMoveToDevice(bonded_particles_lut_manager, 1);

	forcefield = genericMoveToDevice(forcefield, 1);

	cudaDeviceSynchronize();
	printf("Box transferred to device\n");
}

void Box::deleteMembers(const bool is_on_device) {
	if (is_on_device) {
		cudaFree(compounds);
		cudaFree(coordarray_circular_queue);
		cudaFree(solventblockgrid_circular_queue);
		cudaFree(transfermodule_array);

		cudaFree(compound_grid);		
		cudaFree(compound_neighborlists);

		cudaFree(forcefield);

		cudaFree(bridge_bundle);
		cudaFree(bonded_particles_lut_manager);

		/*cudaFree(box->potE_buffer);
		cudaFree(box->traj_buffer);
		cudaFree(box->outdata);
		cudaFree(box->data_GAN);*/
	}
	else {
		delete[] compounds;	// TODO: Finish this
		delete[] bridge_bundle;
		delete[] compound_neighborlists;
		//////delete[] box->solvent_neighborlists;
		delete[] bonded_particles_lut_manager;
	}
}

Box SimUtils::copyToHost(const Box* box_dev) {
	Box box{};
	cudaMemcpy(&box, box_dev, sizeof(Box), cudaMemcpyDeviceToHost);

	//const Compound** cs = &box.compounds;
	genericCopyToHost(&box.compounds, MAX_COMPOUNDS);
	genericCopyToHost(&box.bridge_bundle, 1);

	genericCopyToHost(&box.coordarray_circular_queue, MAX_COMPOUNDS * STEPS_PER_LOGTRANSFER);
	genericCopyToHost(&box.solventblockgrid_circular_queue, STEPS_PER_SOLVENTBLOCKTRANSFER);

	genericCopyToHost(&box.compound_neighborlists, MAX_COMPOUNDS);
	
	genericCopyToHost(&box.bonded_particles_lut_manager, 1);
	genericCopyToHost(&box.forcefield, 1);


	//genericCopyToHost(box.solvent_neighborlists, MAX_COMPOUNDS);
	//box.compounds = genericCopyToHost(box.compounds, box.n_compounds);
	//box.bridge_bundle = genericCopyToHost(box.bridge_bundle, 1);

	//box.coordarray_circular_queue = genericCopyToHost(box.coordarray_circular_queue)
	printf("Box copied to host\n");
	return box;
}

SimulationDevice::SimulationDevice(const SimParams& params_host, std::unique_ptr<Box> box_host) {
	genericCopyToDevice(params_host, &params, 1);
	
	databuffers = new DatabuffersDevice(box_host->total_particles_upperbound, box_host->n_compounds);
	databuffers = genericMoveToDevice(databuffers, 1);

	box_host->moveToDevice();
	cudaMallocManaged(&box, sizeof(Box));
	cudaMemcpy(box, box_host.get(), sizeof(Box), cudaMemcpyHostToDevice);
	box_host.reset();
}

void SimulationDevice::deleteMembers() {
	box->deleteMembers(true);
	cudaFree(box);

	databuffers->freeMembers();
	cudaFree(databuffers);

	cudaFree(params);
}


Simulation::Simulation(const SimParams& ip) :
	simparams_host{ ip }
{
	box_host = std::make_unique<Box>();
}



Simulation::~Simulation() {
	if (sim_dev != nullptr) {
		sim_dev->deleteMembers();
		cudaFree(sim_dev);
	}

}

void Simulation::moveToDevice() {
	if (sim_dev != nullptr) { throw "Expected simdev to be null to move sim to device"; };
	sim_dev = new SimulationDevice(simparams_host, std::move(box_host));
	sim_dev = genericMoveToDevice(sim_dev, 1);
}

void Simulation::copyBoxVariables() {
	n_compounds = box_host->n_compounds;
	n_bridges = box_host->bridge_bundle->n_bridges;


	n_solvents = box_host->n_solvents;
	blocks_per_solventkernel = (int)ceil((float)n_solvents / (float)THREADS_PER_SOLVENTBLOCK);

	compounds_host.resize(n_compounds);
	for (int i = 0; i < n_compounds; i++)
		compounds_host[i] = box_host->compounds[i];

	// Need this variable both on host and device
	total_particles_upperbound = box_host->n_compounds * MAX_COMPOUND_PARTICLES + SolventBlockGrid::blocks_total * MAX_SOLVENTS_IN_BLOCK;
	box_host->total_particles_upperbound = total_particles_upperbound;
}


void InputSimParams::overloadParams(std::map<std::string, double>& dict) {
	overloadParam(dict, &dt, "dt", FEMTO_TO_LIMA);	// convert [fs] to [ls]
	overloadParam(dict, &n_steps, "n_steps");
}

SimParams::SimParams(const InputSimParams& ip) : constparams{ip.n_steps, ip.dt }
{}

DatabuffersDevice::DatabuffersDevice(size_t total_particles_upperbound, int n_compounds) {
	// Permanent Outputs for energy & trajectory analysis
	size_t n_datapoints = total_particles_upperbound * STEPS_PER_LOGTRANSFER;
	printf("Malloc %.2f MB on device for data buffers\n", sizeof(float) * n_datapoints + sizeof(Float3) * n_datapoints);

	cudaMallocManaged(&potE_buffer, sizeof(float) * n_datapoints);
	cudaMallocManaged(&traj_buffer, sizeof(Float3) * n_datapoints);

#ifdef USEDEBUGF3
	uint64_t bytes_for_debugf3 = sizeof(Float3) * DEBUGDATAF3_NVARS * simulation->total_particles_upperbound * simulation->n_steps;
	cudaMallocManaged(&simulation->box->debugdataf3, bytes_for_debugf3);
#endif

	// TRAINING DATA and TEMPRARY OUTPUTS
	int n_loggingdata_device = 10 * STEPS_PER_LOGTRANSFER;
	uint64_t n_traindata_device = static_cast<uint64_t>(N_DATAGAN_VALUES) * MAX_COMPOUND_PARTICLES * n_compounds * STEPS_PER_TRAINDATATRANSFER;
	long double total_bytes = static_cast<long double>(sizeof(float) * static_cast<long double>(n_loggingdata_device) + sizeof(Float3) * n_traindata_device);
	printf("Reserving %.4f MB device mem for logging + training data\n", (float)((total_bytes) * 1e-6));

	cudaMallocManaged(&outdata, sizeof(float) * 10 * STEPS_PER_LOGTRANSFER);	// 10 data streams for 10k steps. 1 step at a time.

	cudaMallocManaged(&data_GAN, sizeof(Float3) * N_DATAGAN_VALUES * MAX_COMPOUND_PARTICLES * n_compounds * STEPS_PER_TRAINDATATRANSFER);
}

void DatabuffersDevice::freeMembers() {
	cudaFree(potE_buffer);
	cudaFree(traj_buffer);

	cudaFree(outdata);
	cudaFree(data_GAN);
}
