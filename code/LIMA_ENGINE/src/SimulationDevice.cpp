#include "SimulationDevice.cuh"
#include "Utilities.h"


//BoxConfig::BoxConfig(Compound* compounds, uint8_t* compoundsAtomTypes, float* compoundsAtomcharges, BondedParticlesLUT* bpLUTs, const BoxGrid::TinymolBlockAdjacency::BlockRef* tinymolNearbyBlockIds, BoxGrid::TinymolBlockAdjacency::NearbyBlocksSequences* tinymolNearbyBlocksSequences) :
//	compounds(compounds),
//	compoundsAtomtypes(compoundsAtomTypes), 
//	compoundsAtomCharges(compoundsAtomcharges),
//	bpLUTs(bpLUTs),
//	tinymolNearbyBlockIds(tinymolNearbyBlockIds),
//	tinymolNearbyBlocksSequences(tinymolNearbyBlocksSequences)
//	//boxparams(boxHost != nullptr ? boxHost->boxparams : BoxParams{}),
//	//uniformElectricField(boxHost != nullptr ? boxHost->uniformElectricField : UniformElectricField{})
//{}
BoxConfig BoxConfig::Create(const Box& boxHost) {
	return BoxConfig();
}
void BoxConfig::FreeMembers() const {

}


BoxState::BoxState(PersistentclusterInterimState* pclusterInterimStates) :
	pclusterInterimStates(pclusterInterimStates)
{}

BoxState BoxState::Create(const Box& boxHost) {
	return BoxState{
		GenericCopyToDevice(boxHost.pclusterInterimStates),
	};
}
void BoxState::CopyDataToHost(Box& boxHost) const {
	cudaMemcpy(boxHost.pclusterInterimStates.data(), pclusterInterimStates, sizeof(PersistentclusterInterimState) * boxHost.pclusterInterimStates.size(), cudaMemcpyDeviceToHost);
}
void BoxState::FreeMembers() const {
	cudaFree(pclusterInterimStates);
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

	potE_buffer = databuffers.potE_buffer;
	traj_buffer = databuffers.traj_buffer;
	vel_buffer = databuffers.vel_buffer;
	forceBuffer = databuffers.forceBuffer;

	/*if (params_host.em_variant) {
		cudaMalloc(&adamState, sizeof(AdamState) * box_host->boxparams.total_particles_upperbound);
		cudaMemset(adamState, 0, sizeof(AdamState) * box_host->boxparams.total_particles_upperbound);
	}*/

	LIMA_UTILS::genericErrorCheck("Error during creation of SimDevice");
}

void SimulationDevice::FreeMembers() {
	boxConfig.FreeMembers();
	boxState.FreeMembers();


	cudaFree(nParticlesInCompoundsBuffer);

	//cudaFree(transfermodule_array);
	cudaFree(signals);


	if (adamState != nullptr)
		cudaFree(adamState);
}

//CompoundQuickData* CompoundQuickData::CreateBuffer(const Simulation& simulation) {
//	std::vector<CompoundQuickData> compoundQuickDataHost(simulation.box_host->boxparams.n_compounds, CompoundQuickData{});
//	memset(compoundQuickDataHost.data(), 0, sizeof(CompoundQuickData) * compoundQuickDataHost.size());
//
//	for (int cid = 0; cid < simulation.box_host->compounds.size(); cid++) {
//		const Compound& compound = simulation.box_host->compounds[cid];
//		CompoundQuickData& quickData = compoundQuickDataHost[cid];
//		for (int pid = 0; pid < MAX_COMPOUND_PARTICLES; pid++) {
//			if (pid < compound.n_particles) {
//				quickData.relPos[pid] = simulation.box_host->compoundCoordsBuffer[cid].rel_positions[pid].ToRelpos();
//				quickData.ljParams[pid] = simulation.forcefield.particle_parameters[compound.atom_types[pid]];
//				quickData.charges[pid] = compound.atom_charges[pid];
//			}
//			else {
//				quickData.relPos[pid] = Float3{};
//				quickData.ljParams[pid] = ForceField_NB::ParticleParameters{};
//				quickData.charges[pid] = 0.f;
//			}
//		}
//	}
//	return GenericCopyToDevice(compoundQuickDataHost);
//}
