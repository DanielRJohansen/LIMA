//#pragma once - this file must NOT be included multiple times

#include "Engine.cuh"
#include "ForceComputations.cuh"
#include "KernelWarnings.cuh"
#include "EngineUtils.cuh"

#include "SimulationDevice.cuh"
#include "BoundaryCondition.cuh"
#include "SolventBlockTransfers.cuh"
#include "DeviceAlgorithms.cuh"

//#include <cuda/pipeline>
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

__global__ void ElasticPositionsForceKernel(const PersistentCluster* const pc, const PersistentClusterMeta* const pcMeta, const Float3* const elasticPositions, ForceEnergy* const forceEnergy, int nPclusters, Float3 boxSize) {
	const int pcId = blockIdx.x * blockDim.x + threadIdx.x;
	if (pcId >= nPclusters)
		return;


	for (int pid = 0; pid < PersistentCluster::maxParticles; pid++) {
		const int pidGlobal = pcMeta[pcId].particleIdsGlobal[pid];
		if (pidGlobal == -1)
			continue;

		const float mass = pcMeta[pcId].mass[pid];
		const Float3 position = pc[pcId].pqd[pid].position;
		const Float3 ep = elasticPositions[pidGlobal];

		Float3 elasticPosition{
			isnan(ep.x) ? position.x : ep.x,
			isnan(ep.y) ? position.y : ep.y,
			isnan(ep.z) ? position.z : ep.z
		};
		PeriodicBoundaryCondition::applyHyperposNM(position, elasticPosition, boxSize);

		const Float3 difference = elasticPosition - position;
		const float distSq = difference.lenSquared();
		const float dist = difference.len();
		const Float3 forceDirection = distSq < 0.00001f ? Float3{ 0.f } : difference.norm();

		
		// Magnitude =  Coeff * (e^(wx^2) - 1) / (e^(wx^2) + 1)  // Coeff controls magnitude, w controls gradient
		const float coefficient = 1000000.f; // [kJ/mol/nm]
		const float w = 6;
		float eTerm = expf(w * dist);
		const float forceMagnitude = mass * coefficient * (eTerm - 1.0f) / (eTerm + 1.0f); 

		const float potentialEnergy = logf(eTerm + 1.0f) * coefficient;

		forceEnergy[pcId * PersistentCluster::maxParticles + pid] = ForceEnergy{ forceDirection * forceMagnitude, potentialEnergy };
	}
}





static const int THREADS_PER_BONDSGROUPSKERNEL = BondGroup::maxParticles;
template <typename BoundaryCondition, bool emVariant>
__global__ void BondgroupsKernel(const BondGroup* const bondGroups, const BoxState boxState, ForceEnergy* const forceEnergiesOut, const PersistentCluster* const pclusters, Float3 boxSize, Float3 boxSizeInv) {
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
		BoundaryCondition::ApplyHyperpos(positions[0], positions[threadIdx.x], boxSize, boxSizeInv);
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
__global__ void NbNonlocalKernel(const SuperCluster* const superClusters, const ScScTask* const tasks, SCResult* const results, const int* const idsOfQuerySuperclusters, const int* const resultIndices, const BoolMatrix16x16* const nointeractionMatrices, 
	const SuperClusterMeta* const superClusterMeta, int step, Float3 boxSize, Float3 boxSizeInv, float ewaldKappa) {
	static_assert(SuperCluster::maxParticles == 16, "This kernel relies on SuperCluster::nParticles being 16");
	__shared__ SuperCluster scSelf;
	__shared__ ScScTask task; // TODO: We dont access this much, no need to store in shared mem...
	//__shared__ ForceEnergy forceenergySelfShared[SuperCluster::maxParticles * 4]; // Blockdim must be 4!
	//__shared__ ForceEnergy feAcc[SuperCluster::maxParticles];

	auto tb = cooperative_groups::this_thread_block();
	cooperative_groups::memcpy_async(tb, &scSelf, &superClusters[blockIdx.x], sizeof(SuperCluster));
	if (threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0) {
		task = tasks[blockIdx.x];
	}	
	//feAcc[threadIdx.x] = ForceEnergy{};
	cooperative_groups::wait(tb);	
	__syncthreads();


	ForceEnergy feInScSelf{};

	const int nBatches = (task.nQueryScs + blockDim.y-1) / blockDim.y;
	for (int batch = 0; batch < nBatches; batch++) {
		const int relativeInteractionIndex = batch * blockDim.y + threadIdx.y;
		const int indexInQueriesBuffer = task.startIndexInQueriesBuffers + relativeInteractionIndex;
		const bool validQuery = relativeInteractionIndex < task.nQueryScs;
		const int queryScId = validQuery ? idsOfQuerySuperclusters[indexInQueriesBuffer] : 0; // For invalid queries we simply load whatever data is at index 0, and continue as normal. This only happens in the final batch, and we dont wanna slow down all other batches with checks

		//cooperative_groups::memcpy_async(tb, &nointeractionsMatrix, &nointeractionMatrices[task.nointeractionMatrixIndex[indexInQueriesBuffer]], sizeof(BoolMatrix16x16));
		PData pdataQueryAtom{};
		superClusters[queryScId].LoadPdata(pdataQueryAtom, threadIdx.x);
		BoundaryCondition::ApplyHyperpos(Float3{ scSelf.posX[0], scSelf.posY[0], scSelf.posZ[0] }, pdataQueryAtom.position, boxSize, boxSizeInv);
		const uint16_t noInteractions = validQuery
			? nointeractionMatrices[indexInQueriesBuffer].GetRow(threadIdx.x)
			: 0xFFFF;

		//ForceEnergy myForceEnergy{};
		ForceEnergy feInQuerySc{};

		for (int i = 0; i < 16; i++) {
			const int indexInScSelf = (threadIdx.x + i) & 15; //% SuperCluster::maxParticles;
			const bool skip = BoolMatrix16x16::Get(noInteractions, indexInScSelf);

			ForceEnergy fe = skip ? ForceEnergy{} : LJ::ComputeParticleParticleNB<computePotE, energyMinimize>(pdataQueryAtom, scSelf, indexInScSelf, -1, -1, ewaldKappa);
			feInQuerySc += fe;

			const int sourceLane = (threadIdx.x - i) & 15;

			fe.force.x = __shfl_sync(0xFFFFFFFFu, fe.force.x, sourceLane, 16);
			fe.force.y = __shfl_sync(0xFFFFFFFFu, fe.force.y, sourceLane, 16);
			fe.force.z = __shfl_sync(0xFFFFFFFFu, fe.force.z, sourceLane, 16);
			fe.potE = __shfl_sync(0xFFFFFFFFu, fe.potE, sourceLane, 16);

			feInScSelf += fe.InvertForce();

			//forceenergySelfShared[threadIdx.y * SuperCluster::maxParticles + indexInScSelf] += fe.InvertForce();
			//__syncwarp();
		}

		if (validQuery) {
			results[resultIndices[indexInQueriesBuffer]].fe[threadIdx.x] = feInQuerySc;
		}
		//__syncthreads(); // Im not sure this is necessary..
	}
	__syncthreads();


	// Reduce self forces
	ForceEnergy* feAcc = (ForceEnergy*)((void*)&scSelf);
	if (threadIdx.y == 0) {
		feAcc[threadIdx.x] = feInScSelf;
	}

	for (int i = 1; i < 4; i++) {
		if (threadIdx.y == i) {
			feAcc[threadIdx.x] += feInScSelf;
		}
		__syncthreads();
	}	
	if (threadIdx.y == 0) {
		const int resultIndex = resultIndices[task.startIndexInQueriesBuffers];
		results[resultIndex].fe[threadIdx.x] = feAcc[threadIdx.x];
	}	
}
 
// blockDim=(16, 4, 1)
template<typename BoundaryCondition, bool emvariant, bool logData>
__global__ void SuperclusterIntegrateKernel(const ForceEnergyInterims forceEnergies, SimulationDevice* const simDev, int data_logging_interval, const SCResult* const scResults,
	SuperCluster* superClusters, const SuperClusterMeta* const scMeta, PersistentCluster* const pclusters, const PersistentClusterMeta* const pcMeta, PersistentclusterInterimState* const pcStates, 
	int64_t step, float dt,	int totalParticlesUpperbound, int numScs, float* forcesMagnitudeSquaredBuffer, /*Only available in EM*/
	Float3 boxSize, float thermostatScalar, Float3* fixedParticleMovementBuffer, Float3* forceMaskBuffer, const Rotation* fixedParticleRotationBuffer /*Only available in LIVEEDIT*/  /*,
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
		p0s[threadIdx.y] = scIdGlobal == -1 ? Float3{} : superClusters[scIdGlobal].Position(0);

		// By applying BC here, we dont need to wait for thread0 later in the kernel
		BoundaryCondition::applyBCNM(p0s[threadIdx.y], boxSize);// TODO: We should use either SC CoM, or a particle close to the middle..
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

	Float3 pos = superClusters[scIdGlobal].Position(threadIdx.x);
	BoundaryCondition::applyHyperposNM(p0s[threadIdx.y], pos, boxSize);

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


	// Write the force to global buffer, only needed to monitor simulations
	forcesMagnitudeSquaredBuffer[pidGlobal] = fe.force.lenSquared();

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
			BoundaryCondition::applyHyperposNM(rotation.center, pos_now, boxSize);
			LAL::RotatePoint(pos_now, rotation.center, rotation.rotation);
			BoundaryCondition::applyHyperposNM(p0s[threadIdx.y], pos_now, boxSize);
			//LAL::RotatePoint(pos_now, Float3{}, Float3{ 0.001f, 0.f, 0.f });
		}
	
		pos = pos_now;// Save pos locally, but only push to box as this kernel ends

		Float3 velScaled;
		velScaled = vel_now * thermostatScalar;

		simDev->boxState.pclusterInterimStates[pcIdGlobal].forces_prev[pidInPcluster] = fe.force;
		simDev->boxState.pclusterInterimStates[pcIdGlobal].vels_prev[pidInPcluster] = velScaled;

		speed = velScaled.len();
	}

	// ------------------------------------------------------------ Boundary Condition --------------------------------------------------------------- //	

	//BoundaryCondition::applyHyperposNM(p0s[threadIdx.y], pos);

	EngineUtils::LogPclusterData(pcIdGlobal, pidInPcluster, step, data_logging_interval, pos, fe.potE, fe.force, speed, totalParticlesUpperbound, simDev);

	superClusters[scIdGlobal].posX[threadIdx.x] = pos.x;
	superClusters[scIdGlobal].posY[threadIdx.x] = pos.y;
	superClusters[scIdGlobal].posZ[threadIdx.x] = pos.z;
	pclusters[pcIdGlobal].pqd[pidInPcluster].position = pos;
}


#pragma warning (pop)
#pragma warning (pop)
