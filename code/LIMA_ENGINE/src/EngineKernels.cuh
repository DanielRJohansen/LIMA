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

		forceEnergy[pcId * PersistentCluster::nParticles + pid] = ForceEnergy{ force, 0.f };
	}
}





static const int THREADS_PER_BONDSGROUPSKERNEL = BondGroup::maxParticles;
template <typename BoundaryCondition, bool emVariant>
__global__ void BondgroupsKernel(const BondGroup* const bondGroups, const BoxState boxState, ForceEnergy* const forceEnergiesOut, const PersistentCluster* const pclusters) {
	__shared__ Float3 positions[BondGroup::maxParticles];

	__shared__ float4 forceEnergyInterrims[BondGroup::maxParticles];

	static const int batchSize = THREADS_PER_BONDSGROUPSKERNEL;
	static const int largestBondBytesize = std::max(sizeof(AngleUreyBradleyBond), sizeof(DihedralBond));
	__shared__ char _bondsBuffer[largestBondBytesize * batchSize];	

	const BondGroup* const bondGroup = &bondGroups[blockIdx.x];
	const BondGroup::ParticleRef pRef = bondGroup->particles[threadIdx.x];	

	forceEnergyInterrims[threadIdx.x] = float4{0,0,0,0};

	// Fetch positions, and hyperpos around first particle.
	if (threadIdx.x < bondGroup->nParticles) {
		positions[threadIdx.x] = pclusters[pRef.pcid].pqd[pRef.pid].position;  //boxState.compoundsRelposNm[pRef.compoundId * MAX_COMPOUND_PARTICLES + pRef.localIdInCompound] + relShift;
	}
	__syncthreads();
	if (threadIdx.x < bondGroup->nParticles) {
		BoundaryCondition::applyHyperposNM(positions[0], positions[threadIdx.x]);
	}
	__syncthreads();

	
	{
		__syncthreads();
		SingleBond* bondsBuffer = reinterpret_cast<SingleBond*>(_bondsBuffer);
		for (int batchStart = 0; batchStart < bondGroup->nSinglebonds; batchStart += blockDim.x) {
			if (batchStart + threadIdx.x < bondGroup->nSinglebonds) {
				const int bondIndex = batchStart + threadIdx.x;
				bondsBuffer[threadIdx.x] = bondGroup->singlebonds[bondIndex];
			}
			__syncthreads();

			LimaForcecalc::computeSinglebondForces<emVariant>(bondsBuffer, std::min(batchSize, bondGroup->nSinglebonds - batchStart), positions, forceEnergyInterrims, 0);
		}
	}

	{
		__syncthreads();
		AngleUreyBradleyBond* bondsBuffer = reinterpret_cast<AngleUreyBradleyBond*>(_bondsBuffer);
		for (int batchStart = 0; batchStart < bondGroup->nAnglebonds; batchStart += blockDim.x) {
			if (batchStart + threadIdx.x < bondGroup->nAnglebonds) {
				const int bondIndex = batchStart + threadIdx.x;
				bondsBuffer[threadIdx.x] = bondGroup->anglebonds[bondIndex];
			}
			__syncthreads();

			LimaForcecalc::computeAnglebondForces(bondsBuffer, std::min(batchSize, bondGroup->nAnglebonds - batchStart), positions, forceEnergyInterrims);
		}
	}

	{
		__syncthreads();
		DihedralBond* bondsBuffer = reinterpret_cast<DihedralBond*>(_bondsBuffer);
		for (int batchStart = 0; batchStart < bondGroup->nDihedralbonds; batchStart += blockDim.x) {
			if (batchStart + threadIdx.x < bondGroup->nDihedralbonds) {
				const int bondIndex = batchStart + threadIdx.x;
				bondsBuffer[threadIdx.x] = bondGroup->dihedralbonds[bondIndex];
			}
			__syncthreads();

			LimaForcecalc::computeDihedralForces(bondsBuffer, std::min(batchSize, bondGroup->nDihedralbonds - batchStart), positions, forceEnergyInterrims);
		}
	}

	{
		__syncthreads();
		ImproperDihedralBond* bondsBuffer = reinterpret_cast<ImproperDihedralBond*>(_bondsBuffer);
		for (int batchStart = 0; batchStart < bondGroup->nImproperdihedralbonds; batchStart += blockDim.x) {
			if (batchStart + threadIdx.x < bondGroup->nImproperdihedralbonds) {
				const int bondIndex = batchStart + threadIdx.x;
				bondsBuffer[threadIdx.x] = bondGroup->improperdihedralbonds[bondIndex];
			}
			__syncthreads();

			LimaForcecalc::computeImproperdihedralForces(bondsBuffer, std::min(batchSize, bondGroup->nImproperdihedralbonds - batchStart), positions, forceEnergyInterrims);
		}
	}


	{
		__syncthreads();
		// TODO: i have no clue if pairbonds should also compute SR electrostatics?
		PairBond* bondsBuffer = reinterpret_cast<PairBond*>(_bondsBuffer);
		for (int batchStart = 0; batchStart < bondGroup->nPairbonds; batchStart += blockDim.x) {
			if (batchStart + threadIdx.x < bondGroup->nPairbonds) {
				const int bondIndex = batchStart + threadIdx.x;
				bondsBuffer[threadIdx.x] = bondGroup->pairbonds[bondIndex];
			}
			__syncthreads();

			LimaForcecalc::computePairbondForces(bondsBuffer, std::min(batchSize, bondGroup->nPairbonds - batchStart), positions, forceEnergyInterrims);
		}
	}

	Float3 force{ forceEnergyInterrims[threadIdx.x].x, forceEnergyInterrims[threadIdx.x].y, forceEnergyInterrims[threadIdx.x].z };
	float potE = forceEnergyInterrims[threadIdx.x].w;

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



// blockdim=16,2,1
template <typename BoundaryCondition, bool energyMinimize, bool computePotE, bool useNointeractionMatrix>
__global__ void NbNonlocalKernel(const SuperCluster* const superClusters, const ScScTask* const tasks, SCResult* const results, const BoolMatrix16x16* const nointeractionMatrices, const SuperClusterMeta* const superClusterMeta, int step) {
	/*__shared__ SuperCluster cluster0;
	__shared__ SuperCluster cluster1;	*/
	__shared__ PData pqd[SuperCluster::nParticles * 2];
	__shared__ ScScTask task;	
	__shared__ Float3 p0Pos; // Used for PBC
	__shared__ BoolMatrix16x16 nointeractionsMatrix;

	const bool controlThread = threadIdx.x == 0 && threadIdx.y == 0;
	uint16_t noInteractions;
	bool hasNoInteractionMatrix = tasks[blockIdx.x].nointeractionMatrixIndex != -1;	

	if (controlThread) {
		task = tasks[blockIdx.x];
		p0Pos = superClusters[task.scIds[0]].pData[0].position;
		if (hasNoInteractionMatrix) {
			nointeractionsMatrix = nointeractionMatrices[task.nointeractionMatrixIndex];
		}
	}
	__syncthreads();

	if (threadIdx.y == 0) {
		pqd[threadIdx.x] = superClusters[task.scIds[0]].pData[threadIdx.x];
		noInteractions = hasNoInteractionMatrix ? nointeractionsMatrix.GetRow(threadIdx.x) : 0;
	}
	else {
		pqd[SuperCluster::nParticles + threadIdx.x] = superClusters[task.scIds[1]].pData[threadIdx.x];
		BoundaryCondition::applyHyperposNM(p0Pos, pqd[SuperCluster::nParticles + threadIdx.x].position); // optim: This reads from __constant__, consider passing the boxSizeHalf directly to the kernel registers??
		noInteractions = hasNoInteractionMatrix ? nointeractionsMatrix.GetColumn(threadIdx.x) : 0;
	}


	ForceEnergy myForceEnergy{};
	const int startQueryIndex = threadIdx.y == 0 ? 16 : 0;
	const int myIndex = threadIdx.y * SuperCluster::nParticles + threadIdx.x;
	for (int i = 0; i < SuperCluster::nParticles; i++) {
		bool skip = false;
		if (hasNoInteractionMatrix && BoolMatrix16x16::Get(noInteractions, i)) {
			skip = true;
		}

		if (!skip) {
			myForceEnergy += LJ::ComputeParticleParticleNB<computePotE, energyMinimize>(pqd[myIndex], pqd[startQueryIndex + i], -1, -1);
		}
	}

	if (threadIdx.y == 0) {
		results[task.resultIndices[0]].fe[threadIdx.x] = myForceEnergy;
	}
	else {
		if (task.scIds[0] != task.scIds[1]) {
			results[task.resultIndices[1]].fe[threadIdx.x] = myForceEnergy;
		}
	}
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
	const int pcIdGlobal = scMeta[blockIdx.x].pclusterIds[threadIdx.x / 4];		// May be -1
	const int pidGlobal = pcIdGlobal == -1 ? -1 : pcMeta[pcIdGlobal].particleIdsGlobal[pidInPcluster];
	positions[threadIdx.x] = superClusters[blockIdx.x].pData[threadIdx.x].position;

	// Collect ForceEnergy from all sources
	ForceEnergy fe{};
	// Gather from NB kernels
	for (int i = scMetaShared.resultsStartIndex; i < scMetaShared.resultsStartIndex + scMetaShared.nResults; i++) {
		KernelHelpersWarnings::ForceCheck(scResults[i].fe[threadIdx.x].force);
		fe += scResults[i].fe[threadIdx.x];	
	}

	fe += pidGlobal == -1 ? ForceEnergy{} : forceEnergies.bonded[pcIdGlobal * PersistentCluster::nParticles + pidInPcluster];
	fe += pidGlobal == -1 ? ForceEnergy{} : forceEnergies.snf[pcIdGlobal * PersistentCluster::nParticles + pidInPcluster];
	fe += pidGlobal == -1 ? ForceEnergy{} : forceEnergies.pme[pcIdGlobal * PersistentCluster::nParticles + pidInPcluster];
	__syncthreads();

	if (pidGlobal != -1) {
		//printf("gid %d bufferIndex %d fx %f\n", pidGlobal, pcIdGlobal * PersistentCluster::nParticles + pidInPcluster, fe.force.x);
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
			const Float3 safeForce = EngineUtils::ForceActivationFunction(fe.force);

			AdamState* const adamState = &simDev->adamState[pcIdGlobal * PersistentCluster::nParticles + pidInPcluster];
			const Float3 pos_now = EngineUtils::IntegratePositionADAM(positions[threadIdx.x], safeForce, adamState, step);
			//printf("posnow %f %f %f\n", pos_now.x, pos_now.y, pos_now.z);

			positions[threadIdx.x] = pos_now;// Save pos locally, but only push to box as this kernel ends
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
