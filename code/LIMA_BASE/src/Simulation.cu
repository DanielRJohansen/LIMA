#include "Simulation.cuh"

Box::~Box() {
	if (owns_members) { deleteMembers(); }
}

void Box::moveToDevice() {
	int bytes_total = sizeof(Compound) * boxparams.n_compounds
		+ sizeof(CompoundState) * MAX_COMPOUNDS * 3
		+ sizeof(NeighborList) * (MAX_SOLVENTS + MAX_COMPOUNDS);
	//printf("BOX: moving %.2f MB to device\n", (float)bytes_total * 1e-6);

	compounds = genericMoveToDevice(compounds, MAX_COMPOUNDS);
	bridge_bundle = genericMoveToDevice(bridge_bundle, 1);

	coordarray_circular_queue = genericMoveToDevice(coordarray_circular_queue, MAX_COMPOUNDS * STEPS_PER_LOGTRANSFER);
	solventblockgrid_circular_queue = genericMoveToDevice(solventblockgrid_circular_queue, STEPS_PER_SOLVENTBLOCKTRANSFER);

	transfermodule_array = genericMoveToDevice(transfermodule_array, SolventBlockGrid::blocks_total);
	compound_grid = genericMoveToDevice(compound_grid, 1);

	compound_neighborlists = genericMoveToDevice(compound_neighborlists, MAX_COMPOUNDS);

	bonded_particles_lut_manager = genericMoveToDevice(bonded_particles_lut_manager, 1);

	forcefield = genericMoveToDevice(forcefield, 1);

	cudaDeviceSynchronize();
	is_on_device = true;
}

void Box::deleteMembers() {
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
		delete[] compounds;
		delete[] coordarray_circular_queue;
		delete[] solventblockgrid_circular_queue;
		delete[] transfermodule_array;

		delete compound_grid;
		delete[] compound_neighborlists;

		delete[] bridge_bundle;
		delete[] bonded_particles_lut_manager;

		// TMP, forcefield should maybe come with other members?
		if (forcefield) { delete forcefield; }

	}
	owns_members = false;	
	cudaDeviceSynchronize();
	auto cuda_status = cudaGetLastError();
	if (cuda_status != cudaSuccess) {
		std::cout << "\nCuda error code: " << cuda_status << " - " << cudaGetErrorString(cuda_status) << std::endl;
		exit(1);
	}
}

std::unique_ptr<Box> SimUtils::copyToHost(Box* box_dev) {
	auto box = std::make_unique<Box>();

	cudaMemcpy(box.get(), box_dev, sizeof(Box), cudaMemcpyDeviceToHost);

	genericCopyToHost(&box->compounds, MAX_COMPOUNDS);
	genericCopyToHost(&box->coordarray_circular_queue, MAX_COMPOUNDS * STEPS_PER_LOGTRANSFER);
	genericCopyToHost(&box->solventblockgrid_circular_queue, STEPS_PER_SOLVENTBLOCKTRANSFER);

	genericCopyToHost(&box->transfermodule_array, SolventBlockGrid::blocks_total);
	genericCopyToHost(&box->compound_grid, 1);

	genericCopyToHost(&box->compound_neighborlists, MAX_COMPOUNDS);
	
	genericCopyToHost(&box->forcefield, 1);

	genericCopyToHost(&box->bridge_bundle, 1);
	genericCopyToHost(&box->bonded_particles_lut_manager, 1);
	

	box->owns_members = true;
	box->is_on_device = false;
	//printf("Box copied to host\n");
	return box;
}

SimulationDevice::SimulationDevice(const SimParams& params_host, std::unique_ptr<Box> box_host) {
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

void SimulationDevice::deleteMembers() {
	box->deleteMembers();
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
	boxparams_host = box_host->boxparams;
	//n_compounds = box_host->boxparams.n_compounds;
	n_bridges = box_host->bridge_bundle->n_bridges;


	//n_solvents = box_host->boxparams.n_solvents;
	//blocks_per_solventkernel = (int)ceil((float)n_solvents / (float)THREADS_PER_SOLVENTBLOCK);

	compounds_host.resize(boxparams_host.n_compounds);
	for (int i = 0; i < boxparams_host.n_compounds; i++)
		compounds_host[i] = box_host->compounds[i];

	// Need this variable both on host and device
	//total_particles_upperbound = box_host->boxparams.n_compounds * MAX_COMPOUND_PARTICLES + SolventBlockGrid::blocks_total * MAX_SOLVENTS_IN_BLOCK;
	//box_host->boxparams.total_particles_upperbound = total_particles_upperbound;
}

void InputSimParams::overloadParams(std::map<std::string, double>& dict) {
	overloadParam(dict, &dt, "dt", FEMTO_TO_LIMA);	// convert [fs] to [ls]
	overloadParam(dict, &n_steps, "n_steps");
}

SimParams::SimParams(const InputSimParams& ip) : constparams{ip.n_steps, ip.dt }
{}

DatabuffersDevice::DatabuffersDevice(size_t total_particles_upperbound, int n_compounds) {
	// Permanent Outputs for energy & trajectory analysis
	{
		const size_t n_datapoints = total_particles_upperbound * STEPS_PER_LOGTRANSFER;
		const size_t bytesize_mb = (sizeof(float) * n_datapoints + sizeof(Float3) * n_datapoints) / 1'000'000;
		assert(n_datapoints && "Tried creating traj or potE buffers with 0 datapoints");
		assert(bytesize_mb < 6'000 && "Tried reserving >6GB data on device");

		cudaMallocManaged(&potE_buffer, sizeof(float) * n_datapoints);
		cudaMallocManaged(&traj_buffer, sizeof(Float3) * n_datapoints);
	}
	

#ifdef USEDEBUGF3
	uint64_t bytes_for_debugf3 = sizeof(Float3) * DEBUGDATAF3_NVARS * simulation->total_particles_upperbound * simulation->n_steps;
	cudaMallocManaged(&simulation->box->debugdataf3, bytes_for_debugf3);
#endif

	// TRAINING DATA and TEMPRARY OUTPUTS
	{
		size_t n_outdata = 10 * STEPS_PER_LOGTRANSFER;
		size_t n_traindata = std::max(static_cast<size_t>(N_DATAGAN_VALUES) * MAX_COMPOUND_PARTICLES * n_compounds * STEPS_PER_TRAINDATATRANSFER, size_t{ 1 });	// This is a hack; can't allocate 0 bytes.
		assert(n_outdata && "Tried to create outdata buffer with 0 datapoints");
		assert(n_traindata && "Tried to create traindata buffer with 0 datapoints");
		
		size_t bytesize_mb = (sizeof(float) * n_outdata + sizeof(Float3) * n_traindata) / 1'000'000;
		assert(bytesize_mb < 6'000 && "Tried reserving >6GB data on device");


		cudaMallocManaged(&outdata, sizeof(float) * n_outdata);	// 10 data streams for 10k steps. 1 step at a time.
		cudaMallocManaged(&data_GAN, sizeof(Float3) * n_traindata);
	}
	
}

void DatabuffersDevice::freeMembers() {
	cudaFree(potE_buffer);
	cudaFree(traj_buffer);

	cudaFree(outdata);
	cudaFree(data_GAN);
}