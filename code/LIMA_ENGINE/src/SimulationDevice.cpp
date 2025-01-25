#include "SimulationDevice.cuh"
#include "Utilities.h"


BoxConfig::BoxConfig(Compound* compounds, uint8_t* compoundsAtomTypes, float* compoundsAtomcharges, BondedParticlesLUT* bpLUTs, const BoxGrid::TinymolBlockAdjacency::BlockRef* tinymolNearbyBlockIds) :
	compounds(compounds),
	compoundsAtomtypes(compoundsAtomTypes), 
	compoundsAtomCharges(compoundsAtomcharges),
	bpLUTs(bpLUTs),
	tinymolNearbyBlockIds(tinymolNearbyBlockIds)
	//boxparams(boxHost != nullptr ? boxHost->boxparams : BoxParams{}),
	//uniformElectricField(boxHost != nullptr ? boxHost->uniformElectricField : UniformElectricField{})
{}
BoxConfig BoxConfig::Create(const Box& boxHost) {
	std::vector<uint8_t> compoundsAtomTypes;
	std::vector<float> compoundsAtomCharges;
	compoundsAtomTypes.reserve(MAX_COMPOUND_PARTICLES * boxHost.boxparams.n_compounds);
	compoundsAtomCharges.reserve(MAX_COMPOUND_PARTICLES * boxHost.boxparams.n_compounds);

	for (int cid = 0; cid < boxHost.boxparams.n_compounds; cid++) { // OPTIM This is very slow
		compoundsAtomTypes.insert(compoundsAtomTypes.end(), boxHost.compounds[cid].atom_types, boxHost.compounds[cid].atom_types + MAX_COMPOUND_PARTICLES);
		compoundsAtomCharges.insert(compoundsAtomCharges.end(), boxHost.compounds[cid].atom_charges, boxHost.compounds[cid].atom_charges + MAX_COMPOUND_PARTICLES);
	}

	return BoxConfig (
		GenericCopyToDevice(boxHost.compounds), 
		GenericCopyToDevice(compoundsAtomTypes),
		GenericCopyToDevice(compoundsAtomCharges), 
		GenericCopyToDevice(boxHost.bpLutCollection), 
        BoxGrid::TinymolBlockAdjacency::PrecomputeNeabyBlockIds(boxHost.boxparams.boxSize, 1.2f)// TODO: MAGIC nr, use the actual cutoff from simparams
	);
}
void BoxConfig::FreeMembers() const {
	BoxConfig boxtemp(nullptr, nullptr, nullptr, nullptr, nullptr);
	cudaMemcpy(&boxtemp, this, sizeof(BoxConfig), cudaMemcpyDeviceToHost);

	cudaFree((void*)boxtemp.compounds);
	cudaFree((void*)boxtemp.compoundsAtomtypes);
	cudaFree((void*)boxtemp.compoundsAtomCharges);
	cudaFree((void*)boxtemp.bpLUTs);
	cudaFree((void*)boxtemp.tinymolNearbyBlockIds);
}


BoxState::BoxState(NodeIndex* compoundOrigos, Float3* compoundsRelpos, TinyMolState* tinyMols,
	SolventBlock* solventblockgrid_circularqueue, CompoundInterimState* compoundsInterimState) :
	compoundOrigos(compoundOrigos), compoundsRelposNm(compoundsRelpos), tinyMols(tinyMols), solventblockgrid_circularqueue(solventblockgrid_circularqueue), compoundsInterimState(compoundsInterimState)
{}
BoxState* BoxState::Create(const Box& boxHost) {
	std::vector<NodeIndex> compoundsOrigos;
	std::vector<Float3> compoundsRelPos;
	for (const auto& compoundCoords : boxHost.compoundCoordsBuffer) {
		compoundsOrigos.push_back(compoundCoords.origo);
		for (int i = 0; i < MAX_COMPOUND_PARTICLES; i++) {
			compoundsRelPos.emplace_back(compoundCoords.rel_positions[i].ToRelpos());
		}
	}

	

	BoxState boxState(		
		GenericCopyToDevice(compoundsOrigos),
		GenericCopyToDevice(compoundsRelPos),
		GenericCopyToDevice(boxHost.tinyMols),
		GenericCopyToDevice(boxHost.solventblockgrid_circularqueue),		
		GenericCopyToDevice(boxHost.compoundInterimStates));

	BoxState* boxStateDev;
	cudaMallocManaged(&boxStateDev, sizeof(BoxState)); // TODO not managed
	cudaMemcpy(boxStateDev, &boxState, sizeof(BoxState), cudaMemcpyHostToDevice);

	return boxStateDev;
}
void BoxState::CopyDataToHost(Box& boxHost) const {
	BoxState boxtemp(nullptr, nullptr, nullptr, nullptr, nullptr);
	cudaMemcpy(&boxtemp, this, sizeof(BoxState), cudaMemcpyDeviceToHost);

	//assert(boxHost.compounds.size() == boxtemp.boxparams.n_compounds);
//	cudaMemcpy(boxHost.compounds.data(), boxtemp.compounds, sizeof(Compound) * boxHost.compounds.size(), cudaMemcpyDeviceToHost); // This should NOT be necessary since the state dont change
	cudaMemcpy(boxHost.compoundInterimStates.data(), boxtemp.compoundsInterimState, sizeof(CompoundInterimState) * boxHost.compoundInterimStates.size(), cudaMemcpyDeviceToHost);
	cudaMemcpy(boxHost.tinyMols.data(), boxtemp.tinyMols, sizeof(TinyMolState) * boxHost.tinyMols.size(), cudaMemcpyDeviceToHost);

	std::vector<NodeIndex> compoundsOrigos;
	std::vector<CompoundInterimState> compoundStates;
	GenericCopyToHost(boxtemp.compoundOrigos, compoundsOrigos, boxHost.compounds.size());
	GenericCopyToHost(boxtemp.compoundsInterimState, compoundStates, boxHost.compounds.size());
	for (int cid = 0; cid < boxHost.compoundCoordsBuffer.size(); cid++) {
		boxHost.compoundCoordsBuffer[cid].origo = compoundsOrigos[cid];
		for (int pid = 0; pid < MAX_COMPOUND_PARTICLES; pid++)
			boxHost.compoundCoordsBuffer[cid].rel_positions[pid] = compoundStates[cid].coords[pid];
	}


	boxHost.solventblockgrid_circularqueue = GenericCopyToHost(boxtemp.solventblockgrid_circularqueue, SolventBlocksCircularQueue::nElementsTotal(boxHost.boxparams.boxSize));	
	LIMA_UTILS::genericErrorCheck("Error during CopyDataToHost\n");
}
void BoxState::FreeMembers() {
	BoxState boxtemp(nullptr, nullptr, nullptr, nullptr, nullptr);
	cudaMemcpy(&boxtemp, this, sizeof(BoxState), cudaMemcpyDeviceToHost);

	cudaFree(boxtemp.compoundsInterimState);
	cudaFree(boxtemp.compoundOrigos);
	cudaFree(boxtemp.compoundsRelposNm);
	cudaFree(boxtemp.tinyMols);
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








SimulationDevice::SimulationDevice(const SimParams& params_host, Box* box_host, const BoxConfig& boxConfig,
	BoxState* boxState, const DatabuffersDeviceController& databuffers) : 
	boxConfig(boxConfig), boxState(boxState), params(params_host),
	boxparams(box_host != nullptr ? box_host->boxparams : BoxParams{})
	//uniformElectricField(box_host != nullptr ? box_host->uniformElectricField : UniformElectricField{})
{
	// Allocate structures for keeping track of solvents and compounds
	compound_grid = BoxGrid::MallocOnDevice<CompoundGridNode>(box_host->boxparams.boxSize);
    cudaMallocManaged(&compound_neighborlists, sizeof(NeighborList) * MAX_COMPOUNDS); // OPTIM hunt down and removed managed mem

	cudaMallocManaged(&transfermodule_array, sizeof(SolventBlockTransfermodule) * BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(box_host->boxparams.boxSize)));

	{
		SimSignals temp{};
		genericCopyToDevice(temp, &signals, 1);
	}

	//genericCopyToDevice(params_host, &params, 1);

    cudaMallocManaged(&nNonbondedNeighborsBuffer, sizeof(uint16_t) * box_host->boxparams.n_compounds);
    cudaMallocManaged(&nonbondedNeighborsBuffer, sizeof(NlistUtil::IdAndRelshift) * box_host->boxparams.n_compounds * NlistUtil::maxCompounds);
    cudaMemset(nNonbondedNeighborsBuffer, 0, sizeof(uint16_t) * box_host->boxparams.n_compounds);


	potE_buffer = databuffers.potE_buffer;
	traj_buffer = databuffers.traj_buffer;
	vel_buffer = databuffers.vel_buffer;
	forceBuffer = databuffers.forceBuffer;

	if (params_host.em_variant) {
		cudaMalloc(&adamState, sizeof(AdamState) * box_host->boxparams.total_particles_upperbound);
		cudaMemset(adamState, 0, sizeof(AdamState) * box_host->boxparams.total_particles_upperbound);
	}

	LIMA_UTILS::genericErrorCheck("Error during creation of SimDevice");
}

void SimulationDevice::FreeMembers() {
	boxConfig.FreeMembers();
	boxState->FreeMembers();
	cudaFree(boxState);



	cudaFree(compound_grid);
	cudaFree(compound_neighborlists);
	cudaFree(transfermodule_array);
	cudaFree(signals);

    cudaFree(nNonbondedNeighborsBuffer);
    cudaFree(nonbondedNeighborsBuffer);

	if (adamState != nullptr)
		cudaFree(adamState);
}

CompoundQuickData* CompoundQuickData::CreateBuffer(const Simulation& simulation) {
	std::vector<CompoundQuickData> compoundQuickDataHost(simulation.box_host->boxparams.n_compounds, CompoundQuickData{});
	for (int cid = 0; cid < simulation.box_host->compounds.size(); cid++) {
		const Compound& compound = simulation.box_host->compounds[cid];
		CompoundQuickData& quickData = compoundQuickDataHost[cid];
		for (int pid = 0; pid < MAX_COMPOUND_PARTICLES; pid++) {
			if (pid < compound.n_particles) {
				quickData.relPos[pid] = simulation.box_host->compoundCoordsBuffer[cid].rel_positions[pid].ToRelpos();
				quickData.ljParams[pid] = simulation.forcefield.particle_parameters[compound.atom_types[pid]];
				quickData.charges[pid] = compound.atom_charges[pid];
			}
			else {
				quickData.relPos[pid] = Float3{};
				quickData.ljParams[pid] = ForceField_NB::ParticleParameters{};
				quickData.charges[pid] = 0.f;
			}
		}
	}
	return GenericCopyToDevice(compoundQuickDataHost);
}
