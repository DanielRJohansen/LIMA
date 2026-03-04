//#pragma once - this file must NOT be included multiple times

#include "Engine.cuh"
#include "ForceComputations.cuh"
#include "KernelWarnings.cuh"
#include "EngineUtils.cuh"

#include "SimulationDevice.cuh"
#include "BoundaryCondition.cuh"
#include "SolventBlockTransfers.cuh"
#include "DeviceAlgorithms.cuh"
#include "Neighborlists.cuh"

//#include <cuda/pipeline>
#include "KernelConstants.cuh"

#include "LennardJonesInteractions.cuh"
#include "ParticleClusters.cuh"

#pragma warning(push)
#pragma warning(disable:E0020)
#pragma warning(push)
#pragma warning(disable: 20054)

#pragma diag_suppress 20054

















// ------------------------------------------------------------------------------------------- KERNELS -------------------------------------------------------------------------------------------//


template <typename BoundaryCondition, bool energyMinimize>
__global__ void PclusterSnfKernel(const PersistentCluster* const pc, const PersistentClusterMeta* const pcMeta, const UniformElectricField uniformElectricField, ForceEnergy* const forceEnergy, int nPclusters) {

	const int pcId = blockIdx.x * blockDim.x + threadIdx.x;	
	if (pcId >= nPclusters)
		return;

	
	for (int pid = 0; pid < PersistentCluster::nParticles; pid++) {
		const int pidGlobal = pcMeta[pcId].particleIdsGlobal[pid];
		if (pidGlobal == -1)
			continue;

		float charge = pc[pcId].pqd[pid].params.charge;
		Float3 force = uniformElectricField.GetForce(charge);

		forceEnergy[pidGlobal] = ForceEnergy{ force, 0.f };
	}
}





static const int THREADS_PER_BONDSGROUPSKERNEL = BondGroup::maxParticles;
template <typename BoundaryCondition, bool emVariant>
__global__ void BondgroupsKernel(const BondGroup* const bondGroups, const BoxState boxState, ForceEnergy* const forceEnergiesOut, const PersistentCluster* const pclusters) {
	__shared__ Float3 positions[BondGroup::maxParticles];

	__shared__ Float3 forcesInterrim[BondGroup::maxParticles];
	__shared__ float potEInterrim[BondGroup::maxParticles];

	static const int batchSize = THREADS_PER_BONDSGROUPSKERNEL;
	static const int largestBondBytesize = std::max(sizeof(AngleUreyBradleyBond), sizeof(DihedralBond));
	__shared__ char _bondsBuffer[largestBondBytesize * batchSize];	

	const BondGroup* const bondGroup = &bondGroups[blockIdx.x];
	const BondGroup::ParticleRef pRef = bondGroup->particles[threadIdx.x];	
	__syncthreads();

	// Fetch positions, and hyperpos around first particle.
	if (threadIdx.x < bondGroup->nParticles) {
		positions[threadIdx.x] = pclusters[pRef.pcid].pqd[pRef.pid].position;  //boxState.compoundsRelposNm[pRef.compoundId * MAX_COMPOUND_PARTICLES + pRef.localIdInCompound] + relShift;
	}
	__syncthreads();
	if (threadIdx.x < bondGroup->nParticles) {
		BoundaryCondition::applyHyperposNM(positions[0], positions[threadIdx.x]);
	}
	__syncthreads();
	

	Float3 force{};
	float potE{};

	
	{
		SingleBond* bondsBuffer = reinterpret_cast<SingleBond*>(_bondsBuffer);
		for (int batchStart = 0; batchStart < bondGroup->nSinglebonds; batchStart += blockDim.x) {
			if (batchStart + threadIdx.x < bondGroup->nSinglebonds) {
				const int bondIndex = batchStart + threadIdx.x;
				bondsBuffer[threadIdx.x] = bondGroup->singlebonds[bondIndex];
			}
			__syncthreads();



			force += LimaForcecalc::computeSinglebondForces<emVariant>(bondsBuffer, std::min(batchSize, bondGroup->nSinglebonds - batchStart), positions, forcesInterrim, potEInterrim, &potE, 0);
		}
	}

	{
		AngleUreyBradleyBond* bondsBuffer = reinterpret_cast<AngleUreyBradleyBond*>(_bondsBuffer);
		for (int batchStart = 0; batchStart < bondGroup->nAnglebonds; batchStart += blockDim.x) {
			if (batchStart + threadIdx.x < bondGroup->nAnglebonds) {
				const int bondIndex = batchStart + threadIdx.x;
				bondsBuffer[threadIdx.x] = bondGroup->anglebonds[bondIndex];
			}
			__syncthreads();

			force += LimaForcecalc::computeAnglebondForces(bondsBuffer, std::min(batchSize, bondGroup->nAnglebonds - batchStart), positions, forcesInterrim, potEInterrim, &potE);
		}
	}

	{
		DihedralBond* bondsBuffer = reinterpret_cast<DihedralBond*>(_bondsBuffer);
		for (int batchStart = 0; batchStart < bondGroup->nDihedralbonds; batchStart += blockDim.x) {
			if (batchStart + threadIdx.x < bondGroup->nDihedralbonds) {
				const int bondIndex = batchStart + threadIdx.x;
				bondsBuffer[threadIdx.x] = bondGroup->dihedralbonds[bondIndex];
			}
			__syncthreads();

			force += LimaForcecalc::computeDihedralForces(bondsBuffer, std::min(batchSize, bondGroup->nDihedralbonds - batchStart), positions, forcesInterrim, potEInterrim, &potE);
		}
	}

	{
		ImproperDihedralBond* bondsBuffer = reinterpret_cast<ImproperDihedralBond*>(_bondsBuffer);
		for (int batchStart = 0; batchStart < bondGroup->nImproperdihedralbonds; batchStart += blockDim.x) {
			if (batchStart + threadIdx.x < bondGroup->nImproperdihedralbonds) {
				const int bondIndex = batchStart + threadIdx.x;
				bondsBuffer[threadIdx.x] = bondGroup->improperdihedralbonds[bondIndex];
			}
			__syncthreads();

			force += LimaForcecalc::computeImproperdihedralForces(bondsBuffer, std::min(batchSize, bondGroup->nImproperdihedralbonds - batchStart), positions, forcesInterrim, potEInterrim, &potE);
		}
	}


	{
		// TODO: i have no clue if pairbonds should also compute SR electrostatics?
		PairBond* bondsBuffer = reinterpret_cast<PairBond*>(_bondsBuffer);
		for (int batchStart = 0; batchStart < bondGroup->nPairbonds; batchStart += blockDim.x) {
			if (batchStart + threadIdx.x < bondGroup->nPairbonds) {
				const int bondIndex = batchStart + threadIdx.x;
				bondsBuffer[threadIdx.x] = bondGroup->pairbonds[bondIndex];
			}
			__syncthreads();

			force += LimaForcecalc::computePairbondForces(bondsBuffer, std::min(batchSize, bondGroup->nPairbonds - batchStart), positions, forcesInterrim, potEInterrim, &potE);
		}
	}

	forceEnergiesOut[blockIdx.x * BondGroup::maxParticles + threadIdx.x] = ForceEnergy{ force, potE };
}

// gridDim = (nPclusters, 1, 1)
// blockDim = (32, 1, 1) // TODO OPTIM: Use y=4, and have 1 particle in pc per y-thread
__global__ void PclusterBondgroupsGather(const PersistentClusterMeta* const pclusterMeta, int nPclusters, const ForceEnergyInterims forceEnergies) {
	const int pcId = blockIdx.x * blockDim.x + threadIdx.x;
	if (pcId >= nPclusters)
		return;

	for (int pid = 0; pid < 4; pid++) {
		int pidGlobal = pclusterMeta[pcId].particleIdsGlobal[pid];
		if (pidGlobal == -1)
			continue;

		BondgroupRefManager beRefs = pclusterMeta[pcId].bondgroupReferences[pid];
		ForceEnergy fe{};
		for (int i = 0; i < beRefs.nBondgroupApperances; i++) {
			BondgroupRef bondgroupRef = beRefs.bondgroupApperances[i];
			fe += forceEnergies.forceEnergiesBondgroups[bondgroupRef.bondgroupId * BondGroup::maxParticles + bondgroupRef.localIndexInBondgroup];
		}

		forceEnergies.bonded[pcId * PersistentCluster::nParticles + pid] = fe;
		//printf("Gatherout pid %d fx %f\n", pidGlobal, fe.force.x);
	}
}

// 
/// <summary>
/// blockdim=16,2,1  OPTIM: This is too complex, for now its just 16,1,1
/// MaskMatrix is either BoolMatrix16x16 or NoMat
/// </summary>
template <typename BoundaryCondition, bool energyMinimize, bool computePotE, bool useNointeractionMatrix>
__global__ void NbNonlocalKernel(const SuperCluster* const superClusters, const ScScTask* const tasks, SCResult* const results, const BoolMatrix16x16* const nointeractionMatrices, const SuperClusterMeta* const superClusterMeta, int step) {
	__shared__ ScScTask task;
	__shared__ SuperCluster queryCluster;
	__shared__ Float3 p0Pos; // Used for PBC`
	__shared__ SCResult utilitySCResult;

	
	uint16_t noInteractionsRow;
	bool hasNoInteractionMatrix = tasks[blockIdx.x].nointeractionMatrixIndex != -1;

	if (threadIdx.x == 0) {
		task = tasks[blockIdx.x];
		queryCluster = superClusters[task.scIds[1]];
	}

	utilitySCResult.fe[threadIdx.x] = ForceEnergy{};
	__syncthreads();
	
	
	if constexpr (useNointeractionMatrix) {
		if (hasNoInteractionMatrix)
			noInteractionsRow = nointeractionMatrices[task.nointeractionMatrixIndex].GetRow(threadIdx.x);
	}
	
	const PData myParticle = superClusters[task.scIds[0]].pData[threadIdx.x];
	if (threadIdx.x == 0) {
		p0Pos = superClusters[task.scIds[0]].pData[0].position;
	}
	__syncthreads();

	BoundaryCondition::applyHyperposNM(p0Pos, queryCluster.pData[threadIdx.x].position); // optim: This reads from __constant__, consider passing the boxSizeHalf directly to the kernel registers??		
	__syncthreads();


	ForceEnergy myForceEnergy{};

	//const int firstQueryIndex = threadIdx.x + threadIdx.y * SuperCluster::nParticles/2;
	for (int i = 0; i < SuperCluster::nParticles; i++) {
		int queryIndex = threadIdx.x + i;
		queryIndex -= SuperCluster::nParticles * (queryIndex >= SuperCluster::nParticles);

		bool skip = false;

		if constexpr (useNointeractionMatrix) {
			if (hasNoInteractionMatrix && BoolMatrix16x16::Get(noInteractionsRow, queryIndex)) {
				skip = true;
			}
		}


		ForceEnergy fe{}; 
		if (!skip) { 
			int p0pid = superClusterMeta[task.scIds[0]].particlesIds[threadIdx.x];
			int p1pid = superClusterMeta[task.scIds[1]].particlesIds[queryIndex];
			if (step != 788) {
				p0pid = -1;
				p1pid = -1;
			}

			fe = LJ::ComputeParticleParticleNB<computePotE, energyMinimize>(myParticle, queryCluster.pData[queryIndex], p0pid, p1pid);
		}
		myForceEnergy += fe;

		// Now invert force and push to query atom also
		fe.force = -fe.force;
		utilitySCResult.fe[queryIndex] += fe;

		__syncthreads();
		/*if (threadIdx.y == 1) {
			queryFE.fe[queryIndex] += fe;
		}
		__syncthreads();*/
	}
	__syncthreads();

	// For selfinteraction tasks all particles have already computed their interactions with queryparticles in this cluster, and thus the queryparticle already has the force from this, meaning we DONT need to push here.
	// I.e. the resultindices are also the same, so we would just be writing the same data again.
	if (task.scIds[0] != task.scIds[1]) {
		auto tb = cooperative_groups::this_thread_block();

		if constexpr (INDEXING_CHECKS){
			if (threadIdx.x == 0 && task.resultIndices[1] < 0)
				printf("Illegal resultindex %d\n", task.resultIndices[1]);
			//printf("forceQuery %f %f %f\n", utilitySCResult.fe[threadIdx.x].force.x, utilitySCResult.fe[threadIdx.x].force.y, utilitySCResult.fe[threadIdx.x].force.z);
		}

		cooperative_groups::memcpy_async(tb, &results[task.resultIndices[1]], &utilitySCResult, sizeof(SCResult));
		cooperative_groups::wait(tb);
	}
	__syncthreads();

	

	// Now reduce across y-dimension in the now vacant queryFE
	utilitySCResult.fe[threadIdx.x] = myForceEnergy;
	//utilitySCResult.fe[threadIdx.x] = ForceEnergy{ Float3{3000.f}, 300.f };
	__syncthreads();
	//if (threadIdx.y == 1) {
	//	utilitySCResult.fe[threadIdx.x] += myForceEnergy;
	//}
	//__syncthreads();

	{
		auto tb = cooperative_groups::this_thread_block();
		cooperative_groups::memcpy_async(tb, &results[task.resultIndices[0]], &utilitySCResult, sizeof(SCResult));
		cooperative_groups::wait(tb);
	}
}


// TODO: This layout can be much smarter
// blockdim = 16,1,1
__global__ void SuperclusterForceenergyReduce(const SuperClusterMeta* const scMeta, const PersistentClusterMeta* const pcMeta, const SCResult* const scResults, ForceEnergy* const particleForceEnergies) {
	__shared__ SuperClusterMeta scMetaShared;
	{
		auto tb = cooperative_groups::this_thread_block();
		cooperative_groups::memcpy_async(tb, &scMetaShared, &scMeta[blockIdx.x], sizeof(SuperClusterMeta));
		cooperative_groups::wait(tb);
	}
	__syncthreads();

	ForceEnergy myFE{};
	int pcId = scMetaShared.pclusterIds[threadIdx.x / 4];
	int pid = threadIdx.x % 4;
	int pidGlobal = pcId == -1 ? -1 : pcMeta[pcId].particleIdsGlobal[threadIdx.x % 4];
	if (pidGlobal == -1)
		return;

	for (int i = scMetaShared.resultsStartIndex; i < scMetaShared.resultsStartIndex + scMetaShared.nResults; i++) {
		/*if (particleId == 4)
			scResults[i].fe[threadIdx.x].force.print('R');	*/

		myFE += scResults[i].fe[threadIdx.x];
		if (isnan(scResults[i].fe[threadIdx.x].force.len()))
			printf("Found nan here %d %d\n", blockIdx.x, threadIdx.x);
	}

	// push

	
	//particleForceEnergies[particleId] = myFE;
	particleForceEnergies[pcId * PersistentCluster::nParticles + pid];
}






 
// blockDim=(16, 1, 1) - 1 warp per supercluster. Todo: use ydimension of 2, so we use 32 threads total
template<typename BoundaryCondition, bool emvariant>
__global__ void SuperclusterIntegrateKernel(const ForceEnergyInterims forceEnergies, SimulationDevice* const simDev, const SCResult* const scResults,
	SuperCluster* superClusters, const SuperClusterMeta* const scMeta, PersistentCluster* const pclusters, const PersistentClusterMeta* pcMeta, PersistentclusterInterimState* const pcStates, int64_t step, float dt,
	int totalParticlesUpperbound) {
	__shared__ Float3 positions[SuperCluster::nParticles];
	__shared__ SuperClusterMeta scMetaShared;

	if (threadIdx.x == 0) {
		scMetaShared = scMeta[blockIdx.x];
	}
	__syncthreads();

	const int pidInPcluster = threadIdx.x % 4;									// Always safe
	//const int pidGlobal = scMeta[blockIdx.x].particlesIds[threadIdx.x];			// May be -1
	const int pcIdGlobal = scMeta[blockIdx.x].pclusterIds[threadIdx.x / 4];		// May be -1
	const int pidGlobal = pcIdGlobal == -1 ? -1 : pcMeta[pcIdGlobal].particleIdsGlobal[pidInPcluster];
	positions[threadIdx.x] = superClusters[blockIdx.x].pData[threadIdx.x].position;

	// Collect ForceEnergy from all sources
	ForceEnergy fe{};
	// Gather from NB kernels
	for (int i = scMetaShared.resultsStartIndex; i < scMetaShared.resultsStartIndex + scMetaShared.nResults; i++) {
		KernelHelpersWarnings::ForceCheck(scResults[i].fe[threadIdx.x].force/*, std::string("SuperclusterIntegrateKernel")*/);
		fe += scResults[i].fe[threadIdx.x];	
	}
	// Gather from bonds : TODO: maybe dont store it ordered like this?
	fe += pidGlobal == -1 ? ForceEnergy{} : forceEnergies.bonded[pcIdGlobal * PersistentCluster::nParticles + pidInPcluster];
	 //TODO: Gather from PME, SNF, others??
	//fe += pidGlobal == -1 ? ForceEnergy{} : forceEnergies.forceEnergySNF[pcIdGlobal * PersistentCluster::nParticles + pidInPcluster];
	__syncthreads();

	if (pidGlobal != -1) {
		//printf("gid %d fx %f\n", pidGlobal, fe.force.x);
		/*fe.force.print('F');
		positions[threadIdx.x].print('P');*/
	}

	// ------------------------------------------------------------ Integration --------------------------------------------------------------- //	
	float speed = 0.f;
	if (pidGlobal != -1) {
		const float mass = pcMeta[pcIdGlobal].mass[pidInPcluster];

		// Energy minimize
		if constexpr (emvariant) {
			// TODO: Handle emvariants/ADAM states
			//const Float3 safeForce = EngineUtils::ForceActivationFunction(fe.force);

			//AdamState* const adamState = &sim->adamState[blockIdx.x * MAX_COMPOUND_PARTICLES + threadIdx.x];
			//const Coord pos_now = EngineUtils::IntegratePositionADAM(compound_coords.rel_positions[threadIdx.x], safeForce, adamState, step);

			//compound_coords.rel_positions[threadIdx.x] = pos_now;// Save pos locally, but only push to box as this kernel ends
		}
		else {			

			const Float3 forcePrev = pcStates[pcIdGlobal].forces_prev[pidInPcluster];
			const Float3 velPrev = pcStates[pcIdGlobal].vels_prev[pidInPcluster];
			const Float3 vel_now = EngineUtils::integrateVelocityVVS(velPrev, forcePrev, fe.force, dt, mass);
			//printf("PC speed %f dt %f force %f mass %f\n", vel_now.len(), dt, fe.force.len(), mass);
			const Float3 pos_now = EngineUtils::IntegratePositionVVS(positions[threadIdx.x], vel_now, fe.force, mass, dt);
			//(pos_now - positions[threadIdx.x]).print('d');
			//pos_now.print('N');
			positions[threadIdx.x] = pos_now;// Save pos locally, but only push to box as this kernel ends
			//compound_coords.rel_positions[threadIdx.x] = pos_now;// Save pos locally, but only push to box as this kernel ends

			Float3 velScaled;
			velScaled = vel_now * DeviceConstants::thermostatScalar;

			simDev->boxState.pclusterInterimStates[pcIdGlobal].forces_prev[pidInPcluster] = fe.force;
			simDev->boxState.pclusterInterimStates[pcIdGlobal].vels_prev[pidInPcluster] = velScaled;

			speed = velScaled.len();
		}
	}
	__syncthreads();

	// ------------------------------------------------------------ Boundary Condition --------------------------------------------------------------- //	
	if (threadIdx.x == 0) {
		BoundaryCondition::applyBCNM(positions[0]);// TODO: We should use either SC CoM, or a particle close to the middle..
	}
	__syncthreads();
	BoundaryCondition::applyHyperposNM(positions[0], positions[threadIdx.x]);

	if (pcIdGlobal != -1) {
		//ParticleToCompoundOrSolventMapping mapping = particleToCompoundOrSolventMapping[pidGlobal];		
		EngineUtils::LogPclusterData(pcIdGlobal, pidInPcluster, step, simDev->params, positions[threadIdx.x], fe.potE, fe.force, speed, totalParticlesUpperbound, simDev);
	}


	// Push positions for next step
	//if (threadIdx.x == 0)
	//	sim->boxState.compoundOrigos[blockIdx.x] = compound_coords.origo;
	//sim->boxState.compoundsInterimState[blockIdx.x].coords[threadIdx.x] = compound_coords.rel_positions[threadIdx.x];
	//sim->boxState.compoundsRelposNm[blockIdx.x * MAX_COMPOUND_PARTICLES + threadIdx.x] = compound_coords.rel_positions[threadIdx.x].ToRelpos();
	//compoundQuickData[blockIdx.x].relPos[threadIdx.x] = compound_coords.rel_positions[threadIdx.x].ToRelpos();
	/*if (pidGlobal != -1) 
		printf("pid %d Pusing pos %f %f %f\n", pidGlobal, positions[threadIdx.x].x, positions[threadIdx.x].y, positions[threadIdx.x].z);*/


	superClusters[blockIdx.x].pData[threadIdx.x].position = positions[threadIdx.x];
	if (pcIdGlobal != -1 )
		pclusters[pcIdGlobal].pqd[pidInPcluster].position = positions[threadIdx.x];
}


#pragma warning (pop)
#pragma warning (pop)
