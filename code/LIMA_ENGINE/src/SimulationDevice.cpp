#include "SimulationDevice.cuh"



BoxConfig::BoxConfig(Compound* compounds, uint8_t* compoundsAtomTypes, half* compoundsAtomcharges, CompoundBridgeBundleCompact* bridgeBundle, BondedParticlesLUT* bpLUTs, const Box* boxHost) :
	compounds(compounds),
	compoundsAtomtypes(compoundsAtomTypes), 
	compoundsAtomCharges(compoundsAtomcharges),
	bridge_bundle(bridgeBundle),
	bpLUTs(bpLUTs),
	boxparams(boxHost != nullptr ? boxHost->boxparams : BoxParams{}),
	uniformElectricField(boxHost != nullptr ? boxHost->uniformElectricField : UniformElectricField{})
{}
BoxConfig* BoxConfig::Create(const Box& boxHost) {
	uint8_t* compoundsAtomtypes;
	half* compoundsAtomCharges;
	cudaMalloc(&compoundsAtomtypes, sizeof(uint8_t) * MAX_COMPOUND_PARTICLES * boxHost.boxparams.n_compounds);
	cudaMalloc(&compoundsAtomCharges, sizeof(half) * MAX_COMPOUND_PARTICLES * boxHost.boxparams.n_compounds);
	for (int cid = 0; cid < boxHost.boxparams.n_compounds; cid++) {
		cudaMemcpy(compoundsAtomtypes + MAX_COMPOUND_PARTICLES * cid, boxHost.compounds[cid].atom_types, sizeof(uint8_t) * MAX_COMPOUND_PARTICLES, cudaMemcpyHostToDevice);
		cudaMemcpy(compoundsAtomCharges + MAX_COMPOUND_PARTICLES * cid, boxHost.compounds[cid].atom_charges, sizeof(half) * MAX_COMPOUND_PARTICLES, cudaMemcpyHostToDevice);
	}

	BoxConfig boxTemp(GenericCopyToDevice(boxHost.compounds), compoundsAtomtypes, compoundsAtomCharges, GenericCopyToDevice(boxHost.bridge_bundle.get(), 1), GenericCopyToDevice(boxHost.bpLutCollection), &boxHost);
	BoxConfig* devPtr;
	cudaMallocManaged(&devPtr, sizeof(BoxConfig));
	cudaMemcpy(devPtr, &boxTemp, sizeof(BoxConfig), cudaMemcpyHostToDevice);

	return devPtr;
}
void BoxConfig::FreeMembers() const {
	BoxConfig boxtemp(nullptr, nullptr, nullptr, nullptr, nullptr, nullptr);
	cudaMemcpy(&boxtemp, this, sizeof(BoxConfig), cudaMemcpyDeviceToHost);

	cudaFree((void*)boxtemp.compoundsAtomtypes);
	cudaFree((void*)boxtemp.compoundsAtomCharges);
	cudaFree((void*)boxtemp.bridge_bundle);
	cudaFree((void*)boxtemp.bpLUTs);
}


BoxState::BoxState(CompoundCoords* compoundcoordsCircularQueue, Solvent* solvents,
	SolventBlocksCircularQueue* solventblockgrid_circularqueue, CompoundInterimState* compoundsInterimState) :
	compoundcoordsCircularQueue(compoundcoordsCircularQueue), solvents(solvents), solventblockgrid_circularqueue(solventblockgrid_circularqueue), compoundsInterimState(compoundsInterimState)
{}
BoxState* BoxState::Create(const Box& boxHost) {
	BoxState boxState(		
		boxHost.compoundcoordsCircularQueue->CopyToDevice(), 
		GenericCopyToDevice(boxHost.solvents), 
		boxHost.solventblockgrid_circularqueue->CopyToDevice(),
		GenericCopyToDevice(boxHost.compoundInterimStates));

	BoxState* boxStateDev;
	cudaMallocManaged(&boxStateDev, sizeof(BoxState)); // TODO not managed
	cudaMemcpy(boxStateDev, &boxState, sizeof(BoxState), cudaMemcpyHostToDevice);

	return boxStateDev;
}
void BoxState::CopyDataToHost(Box& boxHost) const {
	BoxState boxtemp(nullptr, nullptr, nullptr, nullptr);
	cudaMemcpy(&boxtemp, this, sizeof(BoxState), cudaMemcpyDeviceToHost);

	//assert(boxHost.compounds.size() == boxtemp.boxparams.n_compounds);
//	cudaMemcpy(boxHost.compounds.data(), boxtemp.compounds, sizeof(Compound) * boxHost.compounds.size(), cudaMemcpyDeviceToHost); // This should NOT be necessary since the state dont change
	cudaMemcpy(boxHost.compoundInterimStates.data(), boxtemp.compoundsInterimState, sizeof(CompoundInterimState) * boxHost.compoundInterimStates.size(), cudaMemcpyDeviceToHost);
	cudaMemcpy(boxHost.solvents.data(), boxtemp.solvents, sizeof(Solvent) * boxHost.solvents.size(), cudaMemcpyDeviceToHost);

	boxHost.solventblockgrid_circularqueue->CopyDataFromDevice(boxtemp.solventblockgrid_circularqueue);
	boxHost.compoundcoordsCircularQueue->CopyDataFromDevice(boxtemp.compoundcoordsCircularQueue);
}
void BoxState::FreeMembers() {
	BoxState boxtemp(nullptr, nullptr, nullptr, nullptr);
	cudaMemcpy(&boxtemp, this, sizeof(BoxState), cudaMemcpyDeviceToHost);

	cudaFree(boxtemp.compoundcoordsCircularQueue);
	cudaFree(boxtemp.solvents);
	cudaFree(boxtemp.solventblockgrid_circularqueue);
}









DatabuffersDeviceController::DatabuffersDeviceController(int total_particles_upperbound, int n_compounds, int loggingInterval) :
	total_particles_upperbound{ total_particles_upperbound }
{
	// Permanent Outputs for energy & trajectory analysis
	{
		const size_t n_datapoints = total_particles_upperbound * nStepsInBuffer;
		const size_t bytesize_mb = (2 * sizeof(float) * n_datapoints + 2 * sizeof(Float3) * n_datapoints) / 1'000'000;
		assert(n_datapoints && "Tried creating traj or potE buffers with 0 datapoints");
		assert(bytesize_mb < 6'000 && "Tried reserving >6GB data on device");

		cudaMalloc(&potE_buffer, sizeof(*potE_buffer) * n_datapoints);
		cudaMalloc(&traj_buffer, sizeof(*traj_buffer) * n_datapoints);
		cudaMalloc(&vel_buffer, sizeof(*vel_buffer) * n_datapoints);
		cudaMalloc(&forceBuffer, sizeof(*forceBuffer) * n_datapoints);
		

		cudaMemset(potE_buffer, 0, sizeof(float) * n_datapoints);
		cudaMemset(traj_buffer, 0, sizeof(Float3) * n_datapoints);
		cudaMemset(vel_buffer, 0, sizeof(float) * n_datapoints);
		cudaMemset(forceBuffer, 0, sizeof(Float3) * n_datapoints);
	}
}
DatabuffersDeviceController::~DatabuffersDeviceController() {
	cudaFree(potE_buffer);
	cudaFree(traj_buffer);
	cudaFree(vel_buffer);
	cudaFree(forceBuffer);
}








SimulationDevice::SimulationDevice(const SimParams& params_host, Box* box_host, BoxConfig* boxConfig,
	BoxState* boxState, const DatabuffersDeviceController& databuffers) : boxConfig(*boxConfig), boxState(boxState), params(params_host)
{
	// Allocate structures for keeping track of solvents and compounds
	compound_grid = BoxGrid::MallocOnDevice<CompoundGridNode>(box_host->boxparams.boxSize);
	cudaMallocManaged(&compound_neighborlists, sizeof(NeighborList) * MAX_COMPOUNDS);

	cudaMallocManaged(&transfermodule_array, sizeof(SolventBlockTransfermodule) * BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(box_host->boxparams.boxSize)));

	{
		SimSignals temp{};
		genericCopyToDevice(temp, &signals, 1);
	}

	//genericCopyToDevice(params_host, &params, 1);

	potE_buffer = databuffers.potE_buffer;
	traj_buffer = databuffers.traj_buffer;
	vel_buffer = databuffers.vel_buffer;
	forceBuffer = databuffers.forceBuffer;

	chargeGrid = BoxGrid::MallocOnDevice<Electrostatics::ChargeNode>(box_host->boxparams.boxSize);
	chargeGridChargeSums = BoxGrid::MallocOnDevice<float>(box_host->boxparams.boxSize);
	chargeGridOutputForceAndPot = BoxGrid::MallocOnDevice<ForceAndPotential>(box_host->boxparams.boxSize);

	if (params_host.em_variant) {
		cudaMalloc(&adamState, sizeof(AdamState) * box_host->boxparams.total_particles_upperbound);
		cudaMemset(adamState, 0, sizeof(AdamState) * box_host->boxparams.total_particles_upperbound);
	}

	LIMA_UTILS::genericErrorCheck("Error during creation of SimDevice");
}

void SimulationDevice::FreeMembers() {
	boxConfig.FreeMembers();
	//cudaFree((void*)boxConfig);
	boxState->FreeMembers();
	cudaFree(boxState);


	//cudaFree(params);

	cudaFree(compound_grid);
	cudaFree(compound_neighborlists);

	//cudaFree(charge_octtree);
	cudaFree(chargeGrid);
	cudaFree(chargeGridChargeSums);
	cudaFree(chargeGridOutputForceAndPot);

	if (adamState != nullptr)
		cudaFree(adamState);
}

