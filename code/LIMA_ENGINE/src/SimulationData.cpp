#include "SimulationData.h"

#include "Utilities.h"

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
	nParticlesUpperbound{ nPclusters * PersistentCluster::maxParticles }
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








//CompoundQuickData* CompoundQuickData::CreateBuffer(const Simulation& simulation) {
//	std::vector<CompoundQuickData> compoundQuickDataHost(simulation.box->boxparams.n_compounds, CompoundQuickData{});
//	memset(compoundQuickDataHost.data(), 0, sizeof(CompoundQuickData) * compoundQuickDataHost.size());
//
//	for (int cid = 0; cid < simulation.box->compounds.size(); cid++) {
//		const Compound& compound = simulation.box->compounds[cid];
//		CompoundQuickData& quickData = compoundQuickDataHost[cid];
//		for (int pid = 0; pid < MAX_COMPOUND_PARTICLES; pid++) {
//			if (pid < compound.n_particles) {
//				quickData.relPos[pid] = simulation.box->compoundCoordsBuffer[cid].rel_positions[pid].ToRelpos();
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
