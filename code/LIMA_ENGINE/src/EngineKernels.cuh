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

	
	for (int pid = 0; pid < PersistentCluster::maxParticles; pid++) {
		const int pidGlobal = pcMeta[pcId].particleIdsGlobal[pid];
		if (pidGlobal == -1)
			continue;

		float charge = pc[pcId].pqd[pid].params.charge;
		Float3 force = uniformElectricField.GetForce(charge);

		forceEnergy[pcId * PersistentCluster::maxParticles + pid] = ForceEnergy{ force, 0.f };
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

		forceEnergies.bonded[pcId * PersistentCluster::maxParticles + pid] = fe;
		//printf("Gatherout pid %d fx %f\n", pidGlobal, fe.force.x);
	}
}



// blockdim=16,4,1
template <typename BoundaryCondition, bool energyMinimize, bool computePotE, bool useNointeractionMatrix>
__global__ void NbNonlocalKernel(const SuperCluster* const superClusters, const ScScTask* const tasks, SCResult* const results, const BoolMatrix16x16* const nointeractionMatrices, const SuperClusterMeta* const superClusterMeta, int step) {
	static_assert(SuperCluster::maxParticles == 16, "This kernel relies on SuperCluster::nParticles being 16");
	__shared__ PData pqd[SuperCluster::maxParticles * 2];
	__shared__ ScScTask task;	
	__shared__ Float3 p0Pos; // Used for PBC
	__shared__ BoolMatrix16x16 nointeractionsMatrix;
	__shared__ ForceEnergy forceEnergiesShared[SuperCluster::maxParticles * 2];

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
		// Load cluster0
		pqd[threadIdx.x] = superClusters[task.scIds[0]].pData[threadIdx.x];

		// Load cluster1
		pqd[SuperCluster::maxParticles + threadIdx.x] = superClusters[task.scIds[1]].pData[threadIdx.x];
		BoundaryCondition::applyHyperposNM(p0Pos, pqd[SuperCluster::maxParticles + threadIdx.x].position); // optim: This reads from __constant__, consider passing the boxSizeHalf directly to the kernel registers??
	}
	__syncthreads();

	if (threadIdx.y < 2)
		noInteractions = hasNoInteractionMatrix ? nointeractionsMatrix.GetRow(threadIdx.x) : 0;
	else
		noInteractions = hasNoInteractionMatrix ? nointeractionsMatrix.GetColumn(threadIdx.x) : 0;


	ForceEnergy myForceEnergy{};

	const int myIndex = (threadIdx.y > 1) * SuperCluster::maxParticles + threadIdx.x;
	const int startQueryIndex = (SuperCluster::maxParticles + threadIdx.y * SuperCluster::maxParticles/2) % (SuperCluster::maxParticles*2);
	const int queryBase = (threadIdx.y * (SuperCluster::maxParticles / 2)) % SuperCluster::maxParticles;
	for (int i = 0; i < SuperCluster::maxParticles/2; i++) {
		const int queryIndexInCluster = queryBase + i;
		bool skip = false;
		if (hasNoInteractionMatrix && BoolMatrix16x16::Get(noInteractions, queryIndexInCluster)) {
			skip = true;
		}

		if (!skip) {
			myForceEnergy += LJ::ComputeParticleParticleNB<computePotE, energyMinimize>(pqd[myIndex], pqd[startQueryIndex + i], -1, -1);
		}
	}

	if (threadIdx.y == 0)
		forceEnergiesShared[threadIdx.x] = myForceEnergy;
	if (threadIdx.y == 2)
		forceEnergiesShared[threadIdx.x + 16] = myForceEnergy;
	__syncthreads();

	if (threadIdx.y == 1) {
		results[task.resultIndices[0]].fe[threadIdx.x] = forceEnergiesShared[threadIdx.x] + myForceEnergy;
	}
	else if (threadIdx.y == 3) {
		if (task.scIds[0] != task.scIds[1]) {
			results[task.resultIndices[1]].fe[threadIdx.x] = forceEnergiesShared[SuperCluster::maxParticles + threadIdx.x] + myForceEnergy;
		}
	}
}


// blockDim=(16, 4, 1)
//__global__ void NBGather(const SuperClusterMeta* const scMetas, const SCResult* const scResults, ForceEnergy* const feOut /*Sorted by SC, not PC*/) {
//	__shared__ SuperClusterMeta scMetaShared;
//	__shared__ ForceEnergy feShared[SuperCluster::nParticles];
//	__shared__ SCResult scResultsShared[8];
//
//
//	if (threadIdx.x == 0 && threadIdx.y == 0) {
//		scMetaShared = scMetas[blockIdx.x];
//	}
//	__syncthreads();
//
//	ForceEnergy fe{};
//	for (int i = threadIdx.y; i < scMetaShared.nResults; i+=4) {
//		int index = scMetaShared.resultsStartIndex + i;
//		KernelHelpersWarnings::ForceCheck(scResults[index].fe[threadIdx.x].force);
//		fe += scResults[index].fe[threadIdx.x];
//	}
//
//	if (threadIdx.y == 0)
//		feShared[threadIdx.x] = fe;
//	__syncthreads();
//	for (int i = 1; i < 4; i++) {
//		if (threadIdx.y == i)
//			feShared[threadIdx.x] += fe;
//		__syncthreads();
//	}
//	
//	if (threadIdx.y == 0)
//		feOut[blockIdx.x * SuperCluster::nParticles + threadIdx.x] = feShared[threadIdx.x];
//}



 
// blockDim=(16, 4, 1)
template<typename BoundaryCondition, bool emvariant, bool logData>
__global__ void SuperclusterIntegrateKernel(const ForceEnergyInterims forceEnergies, SimulationDevice* const simDev, const SCResult* const scResults,
	SuperCluster* superClusters, const SuperClusterMeta* const scMeta, PersistentCluster* const pclusters, const PersistentClusterMeta* const pcMeta, PersistentclusterInterimState* const pcStates, 
	int64_t step, float dt,	int totalParticlesUpperbound, int numScs, 
	Float3* fixedParticleMovementBuffer, Float3* forceMaskBuffer, const Rotation* fixedParticleRotationBuffer /*Only available in EM*/  /*, 
const ForceEnergy* const nbForceenergy*/) {

	const int nScsPerBlock = 4;

	const int scIdLocal = threadIdx.y;
	const int scIdGlobal = (blockIdx.x * nScsPerBlock + threadIdx.y) < numScs ? (blockIdx.x * nScsPerBlock + threadIdx.y) : -1;
	const int pidLocal = threadIdx.y * SuperCluster::maxParticles + threadIdx.x;


	//__shared__ Float3 positions[SuperCluster::nParticles * nScsPerBlock];
	__shared__ Float3 p0s[nScsPerBlock];
	__shared__ SuperClusterMeta scMetaShared[nScsPerBlock];

	if (threadIdx.x == 0) {
		scMetaShared[threadIdx.y] = scIdGlobal == -1 ? SuperClusterMeta{} : scMeta[scIdGlobal];
		p0s[threadIdx.y] = scIdGlobal == -1 ? Float3{} : superClusters[scIdGlobal].pData[0].position;

		// By applying BC here, we dont need to wait for thread0 later in the kernel
		BoundaryCondition::applyBCNM(p0s[threadIdx.y]);// TODO: We should use either SC CoM, or a particle close to the middle..		
	}
	__syncthreads();

	//const int pidInPcluster = threadIdx.x % 4;									// Always safe
	const int pidInPcluster = scIdGlobal == -1 ? -1 : scMeta[scIdGlobal].indexInPcluster[threadIdx.x];	// May be -1
	const int pcIdGlobal = scIdGlobal == -1 ? -1 :  scMeta[scIdGlobal]._pclusterIds[threadIdx.x];		// May be -1
	const int pidGlobal = scIdGlobal == -1 ? -1 : scMeta[scIdGlobal].globalParticleIds[threadIdx.x];	// May be -1
	//const int pidGlobal = pcIdGlobal == -1 ? -1 : pcMeta[pcIdGlobal].particleIdsGlobal[pidInPcluster];

	if (pidGlobal == -1)
		return;// NO SYNCS AFTER THIS!

	if constexpr (INDEXING_CHECKS) {
		if (pidInPcluster == -1 || pcIdGlobal == -1 || pidGlobal == -1) {
			printf("Indexing check failed in SuperclusterIntegrateKernel! scIdGlobal %d, pidInPcluster %d, pcIdGlobal %d, pidGlobal %d\n", scIdGlobal, pidInPcluster, pcIdGlobal, pidGlobal);
		}
	}

	Float3 pos = superClusters[scIdGlobal].pData[threadIdx.x].position;
	BoundaryCondition::applyHyperposNM(p0s[threadIdx.y], pos);

	// Collect ForceEnergy from all sources
	ForceEnergy fe{};
	// Gather from NB kernels
	for (int i = scMetaShared[scIdLocal].resultsStartIndex; i < scMetaShared[scIdLocal].resultsStartIndex + scMetaShared[scIdLocal].nResults; i++) {
		KernelHelpersWarnings::ForceCheck(scResults[i].fe[threadIdx.x].force);
		fe += scResults[i].fe[threadIdx.x];
	}
//	fe += nbForceenergy[scIdGlobal * SuperCluster::nParticles + threadIdx.x];
	fe += forceEnergies.bonded[pcIdGlobal * PersistentCluster::maxParticles + pidInPcluster];
	fe += forceEnergies.snf[pcIdGlobal * PersistentCluster::maxParticles + pidInPcluster]; // TODO: Should this be compiletime, or maybe launch param decided to ignore? Since most simulations would ignore...
	fe += forceEnergies.pme[pcIdGlobal * PersistentCluster::maxParticles + pidInPcluster]; // TODO: OPTIM: These should follow SC layout, not PC




	// ------------------------------------------------------------ Integration --------------------------------------------------------------- //	
	float speed = 0.f;

	const float mass = pcMeta[pcIdGlobal].mass[pidInPcluster];

	// Energy minimize
	if constexpr (emvariant) {
		const Float3 safeForce = EngineUtils::ForceActivationFunction(fe.force);

		AdamState* const adamState = &simDev->adamState[pcIdGlobal * PersistentCluster::maxParticles + pidInPcluster];
		Float3 pos_now = EngineUtils::IntegratePositionADAM(pos, safeForce, adamState, step);
		//printf("posnow %f %f %f\n", pos_now.x, pos_now.y, pos_now.z);

		// Overrule movement inferred by force, if this value is available AND nonzeory
		/*if (fixedParticleMovementBuffer != nullptr) {
			Float3 fixedMovement = fixedParticleMovementBuffer[pidGlobal];
			if (fixedMovement.lenSquared() > 0) {
				pos_now = pos + fixedMovement;
			}
		}*/

		pos = pos_now;// Save pos locally, but only push to box as this kernel ends
	}
	else {

		if (forceMaskBuffer) {
			fe.force = fe.force * forceMaskBuffer[pidGlobal];
		}

		const Float3 forcePrev = pcStates[pcIdGlobal].forces_prev[pidInPcluster];
		const Float3 velPrev = pcStates[pcIdGlobal].vels_prev[pidInPcluster];
		Float3 vel_now = EngineUtils::integrateVelocityVVS(velPrev, forcePrev, fe.force, dt, mass);
		Float3 pos_now = EngineUtils::IntegratePositionVVS(pos, vel_now, fe.force, mass, dt);

		if constexpr (FORCE_CHECKS) {
			if ((pos_now - pos).len() > 0.5f) {
				printf("Warning: Particle %d in PC %d moved %f nm in one step.\n", pidInPcluster, pcIdGlobal, (pos_now - pos).len());
			}
		}

		// TODO: This should be happening in the EM variant, but i cant get that working properly
		if (fixedParticleMovementBuffer != nullptr) {
			Float3 fixedMovement = fixedParticleMovementBuffer[pidGlobal];
			if (fixedMovement.lenSquared() > 0) {
				fe.force = Float3{};
				vel_now = Float3{};
				pos_now = pos + fixedMovement;
			}				
		}

		if (fixedParticleRotationBuffer != nullptr) {
			fe.force = Float3{};
			vel_now = Float3{}; 			
			const Rotation& rotation = fixedParticleRotationBuffer[pidGlobal];
			BoundaryCondition::applyHyperposNM(rotation.center, pos_now);
			LAL::RotatePoint(pos_now, rotation.center, rotation.rotation);
			BoundaryCondition::applyHyperposNM(p0s[threadIdx.y], pos_now);
			//LAL::RotatePoint(pos_now, Float3{}, Float3{ 0.001f, 0.f, 0.f });
		}
	
		pos = pos_now;// Save pos locally, but only push to box as this kernel ends

		Float3 velScaled;
		velScaled = vel_now * DeviceConstants::thermostatScalar;

		simDev->boxState.pclusterInterimStates[pcIdGlobal].forces_prev[pidInPcluster] = fe.force;
		simDev->boxState.pclusterInterimStates[pcIdGlobal].vels_prev[pidInPcluster] = velScaled;

		speed = velScaled.len();
	}

	// ------------------------------------------------------------ Boundary Condition --------------------------------------------------------------- //	

	//BoundaryCondition::applyHyperposNM(p0s[threadIdx.y], pos);

	EngineUtils::LogPclusterData(pcIdGlobal, pidInPcluster, step, simDev->params, pos, fe.potE, fe.force, speed, totalParticlesUpperbound, simDev);

	superClusters[scIdGlobal].pData[threadIdx.x].position = pos;
	pclusters[pcIdGlobal].pqd[pidInPcluster].position = pos;
}


#pragma warning (pop)
#pragma warning (pop)
