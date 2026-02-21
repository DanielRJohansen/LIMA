#include "SimulationDevice.cuh"
#include "Utilities.h"


BoxConfig::BoxConfig(Compound* compounds, uint8_t* compoundsAtomTypes, float* compoundsAtomcharges, BondedParticlesLUT* bpLUTs, const BoxGrid::TinymolBlockAdjacency::BlockRef* tinymolNearbyBlockIds, BoxGrid::TinymolBlockAdjacency::NearbyBlocksSequences* tinymolNearbyBlocksSequences) :
	compounds(compounds),
	compoundsAtomtypes(compoundsAtomTypes), 
	compoundsAtomCharges(compoundsAtomcharges),
	bpLUTs(bpLUTs),
	tinymolNearbyBlockIds(tinymolNearbyBlockIds),
	tinymolNearbyBlocksSequences(tinymolNearbyBlocksSequences)
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

	return BoxConfig(
		GenericCopyToDevice(boxHost.compounds),
		GenericCopyToDevice(compoundsAtomTypes),
		GenericCopyToDevice(compoundsAtomCharges),
		GenericCopyToDevice(boxHost.bpLutCollection),
		BoxGrid::TinymolBlockAdjacency::PrecomputeNeabyBlockIds(boxHost.boxparams.boxSize, 1.2f),// TODO: MAGIC nr, use the actual cutoff from simparams
		BoxGrid::TinymolBlockAdjacency::PrecomputeNearbyBlockSequences(boxHost.boxparams.boxSize)
	);
}
void BoxConfig::FreeMembers() const {
	//cudaFree((void*)compounds);
	//cudaFree((void*)compoundsAtomtypes);
	//cudaFree((void*)compoundsAtomCharges);
	//cudaFree((void*)bpLUTs);
	//cudaFree((void*)tinymolNearbyBlockIds);
	//cudaFree((void*)tinymolNearbyBlocksSequences);
}


BoxState::BoxState(NodeIndex* compoundsOrigos, Float3* compoundsRelpos, PersistentclusterInterimState* pclusterInterimStates,
	SolventBlock* solventblockgrid_circularqueue, int* nParticlesInSolventblock, int* nParticlesPrefixsumInX
	, ParticleQuickData* solventsParticleQuickdata, ParticleQuickData* solventsParticleQuickDataCompressed
	,BoxGrid::TinymolBlockAdjacency::NearbyBlocksSequencesParticles* tinymolNearbyBlocksSequences
	) :
	compoundOrigos(compoundsOrigos), compoundsRelposNm(compoundsRelpos), pclusterInterimStates(pclusterInterimStates),
	//tinyMolParticlesState(tinyMolParticlesState), 
	solventblockgrid_circularqueue(solventblockgrid_circularqueue), nParticlesInSolventblock(nParticlesInSolventblock), nParticlesPrefixsumInX(nParticlesPrefixsumInX)
	, solventsParticleQuickData(solventsParticleQuickdata), solventsParticleQuickDataCompressed(solventsParticleQuickDataCompressed)
	, tinymolNearbyBlocksSequences(tinymolNearbyBlocksSequences)
	//, solventsRelposNm(solventsRelposNm), solventsAtomtypeIds(solventsAtomtypeIds)
{}

BoxState BoxState::Create(const Box& boxHost) {
	std::vector<NodeIndex> compoundsOrigos;	// OPTIM Initiate with correct size!
	std::vector<Float3> compoundsRelPos;
	//for (const auto& compoundCoords : boxHost.compoundCoordsBuffer) {
	//	compoundsOrigos.emplace_back(compoundCoords.origo);
	//	for (int i = 0; i < MAX_COMPOUND_PARTICLES; i++) {
	//		compoundsRelPos.emplace_back(compoundCoords.rel_positions[i].ToRelpos());
	//	}
	//}

	const size_t nSolventblocks = BoxGrid::BlocksTotal(boxHost.boxparams.boxSize);
	std::vector<int> nParticlesInSolventblock(nSolventblocks, 0);
	/*std::vector<Float3> solventsRelposNm(nSolventblocks * SolventBlock::maxParticles, Float3{});
	std::vector<uint8_t> solventsAtomtypeIds(nSolventblocks * SolventBlock::maxParticles, 0);*/
	std::vector<ParticleQuickData> solventsParticleQuickdata(nSolventblocks * SolventBlock::maxParticles, ParticleQuickData{});
	//for (int i = 0; i < nSolventblocks; i++) {
	//	nParticlesInSolventblock[i] = boxHost.solventblockgrid_circularqueue[i].nParticles;
	//	for (int j = 0; j < boxHost.solventblockgrid_circularqueue[i].nParticles; j++) {
	//		Int3 gridId = BoxGrid::Get3dIndex(i, boxHost.boxparams.boxSize);
	//		solventsParticleQuickdata[i * SolventBlock::maxParticles + j] = ParticleQuickData{
	//			boxHost.solventblockgrid_circularqueue[i].rel_pos[j].ToRelpos(),
	//			{(int8_t)gridId.x, (int8_t)gridId.y, (int8_t)gridId.z},
	//			boxHost.solventblockgrid_circularqueue[i].atomtypeIds[j]
	//		};
	//		/*solventsRelposNm[i * SolventBlock::maxParticles + j] = boxHost.solventblockgrid_circularqueue[i].rel_pos[j].ToRelpos();
	//		solventsAtomtypeIds[i * SolventBlock::maxParticles + j] = boxHost.solventblockgrid_circularqueue[i].atomtypeIds[j];*/
	//	}
	//}

	std::vector<int> nParticlesPrefixsumInX(nSolventblocks, 0);
	{
		for (int z = 0; z < boxHost.boxparams.boxSize.z; z++) {
			for (int y = 0; y < boxHost.boxparams.boxSize.y; y++) {

				int prefixSum = 0;
				for (int x = 0; x < boxHost.boxparams.boxSize.x; x++) {
					int idx = BoxGrid::Get1dIndex(Int3{ x, y, z }, boxHost.boxparams.boxSize);
					nParticlesPrefixsumInX[idx] = prefixSum;
					prefixSum += boxHost.solventblockgrid_circularqueue[idx].nParticles;
				}
			}
		}
	}

	std::vector<BoxGrid::TinymolBlockAdjacency::NearbyBlocksSequencesParticles> tinymolNearbyBlocksSequences(nSolventblocks);

	return BoxState{
		GenericCopyToDevice(compoundsOrigos),
		GenericCopyToDevice(compoundsRelPos),
		GenericCopyToDevice(boxHost.pclusterInterimStates),
		//GenericCopyToDevice(boxHost.tinyMolParticlesState),
		GenericCopyToDevice(boxHost.solventblockgrid_circularqueue),
		//nullptr, nullptr,nullptr
		GenericCopyToDevice(nParticlesInSolventblock),
		GenericCopyToDevice(nParticlesPrefixsumInX),
		GenericCopyToDevice(solventsParticleQuickdata),
		GenericCopyToDevice(solventsParticleQuickdata),
		GenericCopyToDevice(tinymolNearbyBlocksSequences)
		/*GenericCopyToDevice(solventsRelposNm),
		GenericCopyToDevice(solventsAtomtypeIds)*/
	};
}
void BoxState::CopyDataToHost(Box& boxHost) const {
	//BoxState boxtemp( nullptr, nullptr, nullptr, nullptr, nullptr, nullptr);
	//cudaMemcpy(&boxtemp, this, sizeof(BoxState), cudaMemcpyDeviceToHost);

	//assert(boxHost.compounds.size() == boxtemp.boxparams.n_compounds);
//	cudaMemcpy(boxHost.compounds.data(), boxtemp.compounds, sizeof(Compound) * boxHost.compounds.size(), cudaMemcpyDeviceToHost); // This should NOT be necessary since the state dont change
	//cudaMemcpy(boxHost.compoundInterimStates.data(), compoundsInterimState, sizeof(CompoundInterimState) * boxHost.compoundInterimStates.size(), cudaMemcpyDeviceToHost);
	cudaMemcpy(boxHost.pclusterInterimStates.data(), pclusterInterimStates, sizeof(PersistentclusterInterimState) * boxHost.pclusterInterimStates.size(), cudaMemcpyDeviceToHost);
	//cudaMemcpy(boxHost.tinyMolParticlesState.data(), boxtemp.tinyMolParticlesState, sizeof(TinyMolParticleState) * boxHost.tinyMolParticlesState.size(), cudaMemcpyDeviceToHost);

	/*std::vector<NodeIndex> compoundsOrigos;
	std::vector<CompoundInterimState> compoundStates;
	GenericCopyToHost(compoundOrigos, compoundsOrigos, boxHost.compounds.size());
	GenericCopyToHost(compoundsInterimState, compoundStates, boxHost.compounds.size());
	for (int cid = 0; cid < boxHost.compoundCoordsBuffer.size(); cid++) {
		boxHost.compoundCoordsBuffer[cid].origo = compoundsOrigos[cid];
		for (int pid = 0; pid < MAX_COMPOUND_PARTICLES; pid++)
			boxHost.compoundCoordsBuffer[cid].rel_positions[pid] = compoundStates[cid].coords[pid];
	}*/


	//boxHost.solventblockgrid_circularqueue = GenericCopyToHost(solventblockgrid_circularqueue, SolventBlocksCircularQueue::nElementsTotal(boxHost.boxparams.boxSize));	
	//LIMA_UTILS::genericErrorCheck("Error during CopyDataToHost\n");
}
void BoxState::FreeMembers() const {
	//BoxState boxtemp(nullptr, nullptr, nullptr, nullptr, nullptr, nullptr); // TODO No longer necessary, as this is no longer a device ptr
	//cudaMemcpy(&boxtemp, this, sizeof(BoxState), cudaMemcpyDeviceToHost);

	cudaFree(pclusterInterimStates);
	cudaFree(compoundOrigos);
	cudaFree(compoundsRelposNm);
	//cudaFree(boxtemp.tinyMolParticlesState);
	cudaFree(solventblockgrid_circularqueue);

	cudaFree(nParticlesInSolventblock);
	cudaFree(nParticlesPrefixsumInX);
	cudaFree(solventsParticleQuickData);
	cudaFree(solventsParticleQuickDataCompressed);
	cudaFree(tinymolNearbyBlocksSequences);
}









DatabuffersDeviceController::DatabuffersDeviceController(int nPclusters, int loggingInterval) :
	nParticlesUpperbound{ nPclusters * PersistentCluster::nParticles }
{
	// Permanent Outputs for energy & trajectory analysis
	{
		const size_t n_datapoints = nParticlesUpperbound * nStepsInBuffer;
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
	const BoxState& boxState, const DatabuffersDeviceController& databuffers) : 
	boxConfig(boxConfig), boxState(boxState), params(params_host),
	boxparams(box_host != nullptr ? box_host->boxparams : BoxParams{})
	//uniformElectricField(box_host != nullptr ? box_host->uniformElectricField : UniformElectricField{})
{
	//cudaMallocManaged(&transfermodule_array, sizeof(SolventBlockTransfermodule) * BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(box_host->boxparams.boxSize)));

	{
		SimSignals temp{};
		genericCopyToDevice(temp, &signals, 1);
	}

	std::vector<uint8_t> nParticlesInCompoundsVec(boxparams.n_compounds);
	for (int i = 0; i < boxparams.n_compounds; i++) {
		nParticlesInCompoundsVec[i] = box_host->compounds[i].n_particles;
	}
	nParticlesInCompoundsBuffer = GenericCopyToDevice(nParticlesInCompoundsVec);

	std::vector<CompoundInteractionBoundary> compoundInteractionBoundariesVec(boxparams.n_compounds);
	for (int i = 0; i < boxparams.n_compounds; i++) {
		compoundInteractionBoundariesVec[i] = box_host->compounds[i].interaction_boundary;
	}
	compoundsInteractionBoundaryBuffer = GenericCopyToDevice(compoundInteractionBoundariesVec);

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
	boxState.FreeMembers();


	cudaFree(nParticlesInCompoundsBuffer);
	cudaFree(compoundsInteractionBoundaryBuffer);

	//cudaFree(transfermodule_array);
	cudaFree(signals);


	if (adamState != nullptr)
		cudaFree(adamState);
}

CompoundQuickData* CompoundQuickData::CreateBuffer(const Simulation& simulation) {
	std::vector<CompoundQuickData> compoundQuickDataHost(simulation.box_host->boxparams.n_compounds, CompoundQuickData{});
	memset(compoundQuickDataHost.data(), 0, sizeof(CompoundQuickData) * compoundQuickDataHost.size());

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
