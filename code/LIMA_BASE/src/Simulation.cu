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

SimulationDevice::SimulationDevice(const SimParams& params_host, Box* box_host) 
	//: box(box_host)
{
	genericCopyToDevice(params_host, &params, 1);
	
	databuffers = new DatabuffersDevice(box_host->total_particles_upperbound);
	databuffers = genericMoveToDevice(databuffers, 1);

	//box = new Box();
	//box = genericMoveToDevice(box, 1);
	// Move box last so we can use it's host variables
	//box->moveToDevice();
}
void SimulationDevice::deleteMembers() {
	//cudaFree(box);
	cudaFree(databuffers);
	cudaFree(params);
	// Simulation class MUST destroy *this after this destructor has finished
}


Simulation::Simulation(const SimParams& ip) :
	simparams_host{ ip }
{
	box = new Box();
}



Simulation::~Simulation() {
	if (sim_dev != nullptr) {
		sim_dev->deleteMembers();
		cudaFree(sim_dev);
	}

	deleteBoxMembers();	// TODO: move to box destructor
}

void Simulation::moveToDevice() {
	box = genericMoveToDevice(box, 1);
	//genericCopyToDevice(simparams_host, &simparams_device, 1);

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
		cudaFree(box->bridge_bundle);
		cudaFree(box->compound_neighborlists);

		cudaFree(box->bonded_particles_lut_manager);

		/*cudaFree(box->potE_buffer);
		cudaFree(box->traj_buffer);
		cudaFree(box->outdata);
		cudaFree(box->data_GAN);*/

		cudaFree(box);
	}
	else {
		delete[] box->compounds;	// TODO: Finish this
		delete[] box->bridge_bundle;
		delete[] box->compound_neighborlists;
		//////delete[] box->solvent_neighborlists;
		delete[] box->bonded_particles_lut_manager;
		delete box;
	}
}

void InputSimParams::overloadParams(std::map<std::string, double>& dict) {
	overloadParam(dict, &dt, "dt", FEMTO_TO_LIMA);	// convert [fs] to [ls]
	overloadParam(dict, &n_steps, "n_steps");
}

SimParams::SimParams(const InputSimParams& ip) : constparams{ip.n_steps, ip.dt }
{}

DatabuffersDevice::DatabuffersDevice(size_t total_particles_upperbound) {
	// Permanent Outputs for energy & trajectory analysis
	size_t n_datapoints = total_particles_upperbound * STEPS_PER_LOGTRANSFER;
	printf("Malloc %.2f MB on device for data buffers\n", sizeof(float) * n_datapoints + sizeof(Float3) * n_datapoints);

	cudaMallocManaged(&potE_buffer, sizeof(float) * n_datapoints);
	cudaMallocManaged(&traj_buffer, sizeof(Float3) * n_datapoints);
}

//DatabuffersDevice::~DatabuffersDevice() {
//	cudaFree(potE_buffer);
//	cudaFree(traj_buffer);
//}