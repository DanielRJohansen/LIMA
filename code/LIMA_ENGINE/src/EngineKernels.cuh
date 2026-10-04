//#pragma once - this file must NOT be included multiple times

#include "Engine.cuh"
#include "ForceComputations.cuh"
#include "KernelWarnings.cuh"
#include "EngineUtils.cuh"

#include "BoundaryCondition.cuh"
#include "DeviceAlgorithms.cuh"
#include "Utilities.h"

//#include <cuda/pipeline>
#include "LennardJonesInteractions.cuh"



















// ------------------------------------------------------------------------------------------- KERNELS -------------------------------------------------------------------------------------------//


template <typename BoundaryCondition, bool energyMinimize>
__global__ void PclusterSnfKernel(const PersistentCluster* const pc, const PersistentClusterMeta* const pcMeta, UniformElectricField electricField, ForceEnergy* const forceEnergy, int pclusterOffset, int nPclusters) {

	const int workId = blockIdx.x * blockDim.x + threadIdx.x;
	if (workId >= nPclusters) return;
	const int pcId = pclusterOffset + workId;

	
	for (int pid = 0; pid < PersistentCluster::maxParticles; pid++) {
		const int pidGlobal = pcMeta[pcId].particleIdsGlobal[pid];
		if (pidGlobal == -1)
			continue;

		float charge = pc[pcId].pqd[pid].params.charge;
		Float3 force = electricField.GetForce(charge);

		forceEnergy[pcId * PersistentCluster::maxParticles + pid] = ForceEnergy{ force, 0.f };
	}
}

__global__ void ElasticPositionsForceKernel(const PersistentCluster* const pc, const PersistentClusterMeta* const pcMeta, const Float3* const elasticPositions, ForceEnergy* const forceEnergy, int pclusterOffset, int nPclusters, Float3 boxSize) {
	const int workId = blockIdx.x * blockDim.x + threadIdx.x;
	if (workId >= nPclusters) return;
	const int pcId = pclusterOffset + workId;


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





static const int THREADS_PER_BONDSGROUPSKERNEL = 64;
template <typename BoundaryCondition, bool emVariant>
__global__ void BondgroupsKernel(const BondGroupsDevice bondGroups, const BoxState boxState, ForceEnergy* const forceEnergiesOut, const PersistentCluster* const pclusters, Float3 boxSize, Float3 boxSizeInv) {
	__shared__ Float3 positions[THREADS_PER_BONDSGROUPSKERNEL];

	__shared__ float4 forceEnergyInterrims[THREADS_PER_BONDSGROUPSKERNEL];

	static const int batchSize = THREADS_PER_BONDSGROUPSKERNEL;
	static const int largestBondBytesize = std::max(sizeof(AngleUreyBradleyBond), sizeof(DihedralBond));
	__shared__ char _bondsBuffer[largestBondBytesize * batchSize];	

	const BondGroup* const bondGroup = &bondGroups.groups[blockIdx.x];

	forceEnergyInterrims[threadIdx.x] = float4{0,0,0,0};

	// Fetch positions. Periodic boundaries are applied per bond, since a group may contain several molecules far apart
	if (threadIdx.x < bondGroup->nParticles) {
		const BondGroup::ParticleRef pRef = bondGroups.particles[bondGroup->indexOfFirstParticle + threadIdx.x];
		positions[threadIdx.x] = pclusters[pRef.pcid].pqd[pRef.pid].position;  //boxState.compoundsRelposNm[pRef.compoundId * MAX_COMPOUND_PARTICLES + pRef.localIdInCompound] + relShift;
	}
	__syncthreads();

	
	{
		__syncthreads();
		SingleBond* bondsBuffer = reinterpret_cast<SingleBond*>(_bondsBuffer);
		for (int batchStart = 0; batchStart < bondGroup->nSinglebonds; batchStart += blockDim.x) {
			if (batchStart + threadIdx.x < bondGroup->nSinglebonds) {
				const int bondIndex = batchStart + threadIdx.x;
				bondsBuffer[threadIdx.x] = bondGroups.singlebonds[bondGroup->indexOfFirstSinglebond + bondIndex];
			}
			__syncthreads();

			LimaForcecalc::computeSinglebondForces<BoundaryCondition, emVariant>(bondsBuffer, std::min(batchSize, bondGroup->nSinglebonds - batchStart), positions, forceEnergyInterrims, 0, boxSize, boxSizeInv);
		}
	}

	{
		__syncthreads();
		AngleUreyBradleyBond* bondsBuffer = reinterpret_cast<AngleUreyBradleyBond*>(_bondsBuffer);
		for (int batchStart = 0; batchStart < bondGroup->nAnglebonds; batchStart += blockDim.x) {
			if (batchStart + threadIdx.x < bondGroup->nAnglebonds) {
				const int bondIndex = batchStart + threadIdx.x;
				bondsBuffer[threadIdx.x] = bondGroups.anglebonds[bondGroup->indexOfFirstAnglebond + bondIndex];
			}
			__syncthreads();

			LimaForcecalc::computeAnglebondForces<BoundaryCondition>(bondsBuffer, std::min(batchSize, bondGroup->nAnglebonds - batchStart), positions, forceEnergyInterrims, boxSize, boxSizeInv);
		}
	}

	{
		__syncthreads();
		DihedralBond* bondsBuffer = reinterpret_cast<DihedralBond*>(_bondsBuffer);
		for (int batchStart = 0; batchStart < bondGroup->nDihedralbonds; batchStart += blockDim.x) {
			if (batchStart + threadIdx.x < bondGroup->nDihedralbonds) {
				const int bondIndex = batchStart + threadIdx.x;
				bondsBuffer[threadIdx.x] = bondGroups.dihedralbonds[bondGroup->indexOfFirstDihedralbond + bondIndex];
			}
			__syncthreads();

			LimaForcecalc::computeDihedralForces<BoundaryCondition>(bondsBuffer, std::min(batchSize, bondGroup->nDihedralbonds - batchStart), positions, forceEnergyInterrims, boxSize, boxSizeInv);
		}
	}

	{
		__syncthreads();
		ImproperDihedralBond* bondsBuffer = reinterpret_cast<ImproperDihedralBond*>(_bondsBuffer);
		for (int batchStart = 0; batchStart < bondGroup->nImproperdihedralbonds; batchStart += blockDim.x) {
			if (batchStart + threadIdx.x < bondGroup->nImproperdihedralbonds) {
				const int bondIndex = batchStart + threadIdx.x;
				bondsBuffer[threadIdx.x] = bondGroups.improperdihedralbonds[bondGroup->indexOfFirstImproperdihedralbond + bondIndex];
			}
			__syncthreads();

			LimaForcecalc::computeImproperdihedralForces<BoundaryCondition>(bondsBuffer, std::min(batchSize, bondGroup->nImproperdihedralbonds - batchStart), positions, forceEnergyInterrims, boxSize, boxSizeInv);
		}
	}


	{
		__syncthreads();
		// TODO: i have no clue if pairbonds should also compute SR electrostatics?
		PairBond* bondsBuffer = reinterpret_cast<PairBond*>(_bondsBuffer);
		for (int batchStart = 0; batchStart < bondGroup->nPairbonds; batchStart += blockDim.x) {
			if (batchStart + threadIdx.x < bondGroup->nPairbonds) {
				const int bondIndex = batchStart + threadIdx.x;
				bondsBuffer[threadIdx.x] = bondGroups.pairbonds[bondGroup->indexOfFirstPairbond + bondIndex];
			}
			__syncthreads();

			LimaForcecalc::computePairbondForces<BoundaryCondition>(bondsBuffer, std::min(batchSize, bondGroup->nPairbonds - batchStart), positions, forceEnergyInterrims, boxSize, boxSizeInv);
		}
	}

	Float3 force{ forceEnergyInterrims[threadIdx.x].x, forceEnergyInterrims[threadIdx.x].y, forceEnergyInterrims[threadIdx.x].z };
	float potE = forceEnergyInterrims[threadIdx.x].w;

	if (threadIdx.x < bondGroup->nParticles)
		forceEnergiesOut[bondGroup->indexOfFirstParticle + threadIdx.x] = ForceEnergy{ force, potE };
}

// gridDim = (nPclusters, 1, 1)
// blockDim = (32, 1, 1) // TODO OPTIM: Use y=4, and have 1 particle in pc per y-thread
__global__ void PclusterBondgroupsGather(const PersistentClusterMeta* const pclusterMeta, int nPclusters, const ForceEnergyInterims forceEnergies) {
	const int workId = blockIdx.x * blockDim.x + threadIdx.x;
	if (workId >= nPclusters) return;
	const int pcId = workId;

	for (int pid = 0; pid < 4; pid++) {
		int pidGlobal = pclusterMeta[pcId].particleIdsGlobal[pid];
		if (pidGlobal == -1)
			continue;

		BondgroupRefManager beRefs = pclusterMeta[pcId].bondgroupReferences[pid];
		ForceEnergy fe{};
		for (int i = 0; i < beRefs.nBondgroupApperances; i++) {
			BondgroupRef bondgroupRef = beRefs.bondgroupApperances[i];
			fe += forceEnergies.forceEnergiesBondgroups[bondgroupRef.indexInForceEnergiesBondgroups];
		}

		forceEnergies.bonded[pcId * PersistentCluster::maxParticles + pid] = fe;
		//printf("Gatherout pid %d fx %f\n", pidGlobal, fe.force.x);
	}
}



// Deterministic accumulation of NB forces, replacing SCResult for MD: each partial force is converted to 64-bit fixed point
// and summed with integer atomics. Unlike float atomics, integer addition is associative, so the sum is bitwise
// independent of the order the atomics arrive in. The accumulator is small enough to stay L2 resident.
// Not used in EM, where forces can exceed the fixed point range.
struct NbForceAccumulator {
	static constexpr float scale = 16777216.f;		// 2^24 -> resolution 6e-8 J/mol/nm, range +-5.5e11 J/mol/nm
	static constexpr float scaleInv = 1.f / scale;	// Power of 2, so ToFloat is exactly the rounded sum

	// SoA, indexed by scId * SuperCluster::maxParticles + particleIndex. potE is only zeroed/used on logging steps
	unsigned long long* fx = nullptr;
	unsigned long long* fy = nullptr;
	unsigned long long* fz = nullptr;
	unsigned long long* potE = nullptr;

	__device__ static unsigned long long ToFixed(float v) { return static_cast<unsigned long long>(llrintf(v * scale)); }
	__device__ static float ToFloat(unsigned long long v) { return static_cast<float>(static_cast<long long>(v)) * scaleInv; }

	template <bool withPotE>
	__device__ void Add(int index, const ForceEnergy& fe) const {
		atomicAdd(&fx[index], ToFixed(fe.force.x));
		atomicAdd(&fy[index], ToFixed(fe.force.y));
		atomicAdd(&fz[index], ToFixed(fe.force.z));
		if constexpr (withPotE)
			atomicAdd(&potE[index], ToFixed(fe.potE));
	}

	template <bool withPotE>
	__device__ ForceEnergy Load(int index) const {
		return ForceEnergy{ Float3{ ToFloat(fx[index]), ToFloat(fy[index]), ToFloat(fz[index]) }, withPotE ? ToFloat(potE[index]) : 0.f };
	}
};

// blockdim=16,4,1
// computePotE must match logData of the following SuperclusterIntegrateKernel, as results only contain potE when computePotE
//
// __launch_bounds__(64, 20): max 64 threads per block, and we want at least 20 blocks resident per SM.
// This kernel is bound by instruction issue, so it needs many resident warps to always have one ready to issue.
// Residency is limited by the SM's 64K register file: 65536 / (20 blocks * 64 threads) = 51 -> 48 registers per thread.
// Without the bound the compiler happily picks ~59 registers (rounded to 64), which only fits 16 blocks and was measurably slower.
// - The kernel must never be launched with more than 64 threads per block, or the launch fails
// - If future changes need more than 48 registers, the compiler spills to (slow) local memory instead of growing,
//   so after larger changes check ncu for local memory traffic, and revisit the bound
template <typename BoundaryCondition, bool energyMinimize, bool computePotE, bool useNointeractionMatrix>
__global__ void __launch_bounds__(64, 20) NbNonlocalKernel(const SuperCluster* const superClusters, const ScScTask* const tasks,
	SCResult* const results /*Only used in EM*/, const NbForceAccumulator nbForceAcc /*Only used in MD*/,
	const int* const idsOfQuerySuperclusters, const int* const resultIndices, const BoolMatrix16x16* const nointeractionMatrices,
	const SuperClusterMeta* const superClusterMeta, Float3 boxSize, Float3 boxSizeInv, float ewaldKappa) {
	const int scId = blockIdx.x;
	static_assert(SuperCluster::maxParticles == 16, "This kernel relies on SuperCluster::nParticles being 16");
	__shared__ SuperCluster scSelf;
	__shared__ ScScTask task; // TODO: We dont access this much, no need to store in shared mem...
	//__shared__ ForceEnergy forceenergySelfShared[SuperCluster::maxParticles * 4]; // Blockdim must be 4!
	//__shared__ ForceEnergy feAcc[SuperCluster::maxParticles];

	auto tb = cooperative_groups::this_thread_block();
	cooperative_groups::memcpy_async(tb, &scSelf, &superClusters[scId], sizeof(SuperCluster));
	if (threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0) {
		task = tasks[scId];
	}	
	//feAcc[threadIdx.x] = ForceEnergy{};
	cooperative_groups::wait(tb);
	__syncthreads();

	if constexpr (!energyMinimize) {
		if (threadIdx.y == 0)
			LJ::PrescaleNBParams(scSelf.epsilonSqrt[threadIdx.x], scSelf.charge[threadIdx.x]);
		__syncthreads();
	}


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
		if constexpr (!energyMinimize)
			LJ::PrescaleNBParams(pdataQueryAtom.params.epsilonSqrt, pdataQueryAtom.params.charge);
		BoundaryCondition::ApplyHyperpos(Float3{ scSelf.posX[0], scSelf.posY[0], scSelf.posZ[0] }, pdataQueryAtom.position, boxSize, boxSizeInv);
		const uint16_t noInteractions = validQuery
			? nointeractionMatrices[indexInQueriesBuffer].GetRow(threadIdx.x)
			: 0xFFFF;

		//ForceEnergy myForceEnergy{};
		ForceEnergy feInQuerySc{};

		for (int i = 0; i < 16; i++) {
			const int indexInScSelf = (threadIdx.x + i) & 15; //% SuperCluster::maxParticles;
			const bool masked = BoolMatrix16x16::Get(noInteractions, indexInScSelf);

			ForceEnergy fe = LJ::ComputeParticleParticleNB<computePotE, energyMinimize>(pdataQueryAtom, scSelf, indexInScSelf, masked, ewaldKappa);
			feInQuerySc += fe;

			const int sourceLane = (threadIdx.x - i) & 15;

			fe.force.x = __shfl_sync(0xFFFFFFFFu, fe.force.x, sourceLane, 16);
			fe.force.y = __shfl_sync(0xFFFFFFFFu, fe.force.y, sourceLane, 16);
			fe.force.z = __shfl_sync(0xFFFFFFFFu, fe.force.z, sourceLane, 16);
			if constexpr (computePotE)
				fe.potE = __shfl_sync(0xFFFFFFFFu, fe.potE, sourceLane, 16);

			feInScSelf += fe.InvertForce();

			//forceenergySelfShared[threadIdx.y * SuperCluster::maxParticles + indexInScSelf] += fe.InvertForce();
			//__syncwarp();
		}

		if (validQuery) {
			if constexpr (energyMinimize)
				results[resultIndices[indexInQueriesBuffer]].Store<computePotE>(threadIdx.x, feInQuerySc);
			// The self-task's query forces are dropped: its pairs are evaluated in both orders, so the reaction forces in feInScSelf
			// already hold the full force. (In the EM path this result is overwritten by the self reduction below)
			else if (queryScId != scId)
				nbForceAcc.Add<computePotE>(queryScId * SuperCluster::maxParticles + threadIdx.x, feInQuerySc);
		}
		//__syncthreads(); // Im not sure this is necessary..
	}

	if constexpr (!energyMinimize) {
		// No need to reduce the rows first, the fixed point sum is order independent
		nbForceAcc.Add<computePotE>(scId * SuperCluster::maxParticles + threadIdx.x, feInScSelf);
		return;
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
		results[resultIndex].Store<computePotE>(threadIdx.x, feAcc[threadIdx.x]);
	}	
}
 
// blockDim=(16, 4, 1)
template<typename BoundaryCondition, bool emvariant, bool logData>
__global__ void SuperclusterIntegrateKernel(const ForceEnergyInterims forceEnergies, Float3* const emForces /*Only available in EM*/, int data_logging_interval, 
	const SCResult* const scResults /*Only used in EM*/, const NbForceAccumulator nbForceAcc /*Only used in MD*/,
	SuperCluster* superClusters, const SuperClusterMeta* const scMeta, PersistentCluster* const pclusters, const PersistentClusterMeta* const pcMeta, PersistentclusterInterimState* const pcStates, 
	int64_t step, const IntegrationSimulationData* simulationData, int nSuperclusters, float* forcesMagnitudeSquaredBuffer, /*Only available in EM*/
	Float3 boxSize, Float3* fixedParticleMovementBuffer, Float3* forceMaskBuffer, const Rotation* fixedParticleRotationBuffer,
	Float3* trajBuffer, float* potEBuffer, float* velocityBuffer, Float3* forceBuffer /*Only available in LIVEEDIT*/  /*,
const ForceEnergy* const nbForceenergy*/) {

	const int nScsPerBlock = 4;

	const int scIdLocal = threadIdx.y;
	const int scIdGlobal = (blockIdx.x * nScsPerBlock + threadIdx.y) < nSuperclusters ? static_cast<int>(blockIdx.x * nScsPerBlock + threadIdx.y) : -1;


	//__shared__ Float3 positions[SuperCluster::nParticles * nScsPerBlock];
	__shared__ Float3 p0s[nScsPerBlock];
	__shared__ int resultsStartIndex[nScsPerBlock];
	__shared__ int nResults[nScsPerBlock];
	__shared__ IntegrationSimulationData simulations[nScsPerBlock];


	if (threadIdx.x == 0) {
		resultsStartIndex[threadIdx.y] = scIdGlobal == -1 ? -1 : scMeta[scIdGlobal].resultsStartIndex;
		nResults[threadIdx.y] = scIdGlobal == -1 ? 0 : scMeta[scIdGlobal].nResults;
		simulations[threadIdx.y] = scIdGlobal == -1 ? IntegrationSimulationData{} : simulationData[scMeta[scIdGlobal].simulationId];
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
	const auto& simulation = simulations[scIdLocal];

	// Collect ForceEnergy from all sources
	ForceEnergy fe{};
	// Gather from NB kernels
	if constexpr (emvariant) {
		for (int i = resultsStartIndex[scIdLocal]; i < resultsStartIndex[scIdLocal] + nResults[scIdLocal]; i++) {
			const ForceEnergy result = scResults[i].Load<logData>(threadIdx.x);
			KernelHelpersWarnings::ForceCheck(result.force);
			fe += result;
		}
	}
	else {
		const ForceEnergy result = nbForceAcc.Load<logData>(scIdGlobal * SuperCluster::maxParticles + threadIdx.x);
		KernelHelpersWarnings::ForceCheck(result.force);
		fe += result;
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

	// Energy minimize. The particles are moved by EM::UpdateKernel once the step is decided for the whole simulation
	if constexpr (emvariant) {
		// Overlapping particles in unminimized structures can produce infinite forces. Push them apart in a pseudorandom
		// direction, with a force large enough to still dominate the max force
		const bool finite = isfinite(fe.force.lenSquared());
		emForces[pcIdGlobal * PersistentCluster::maxParticles + pidInPcluster] = finite ? fe.force : EngineUtils::GenerateRandomForce(pidGlobal) * 1e9f;
	}
	else {

		if (forceMaskBuffer) {
			fe.force = fe.force * forceMaskBuffer[pidGlobal];
		}

		const Float3 forcePrev = pcStates[pcIdGlobal].forces_prev[pidInPcluster];
		const Float3 velPrev = pcStates[pcIdGlobal].vels_prev[pidInPcluster];
		Float3 vel_now = EngineUtils::integrateVelocityVVS(velPrev, forcePrev, fe.force, simulation.dt, mass);
		Float3 pos_now = EngineUtils::IntegratePositionVVS(pos, vel_now, fe.force, mass, simulation.dt);

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
		velScaled = vel_now * simulation.thermostatScalar;

		pcStates[pcIdGlobal].forces_prev[pidInPcluster] = fe.force;
		pcStates[pcIdGlobal].vels_prev[pidInPcluster] = velScaled;

		speed = velScaled.len();
	}

	// ------------------------------------------------------------ Boundary Condition --------------------------------------------------------------- //	

	//BoundaryCondition::applyHyperposNM(p0s[threadIdx.y], pos);

	EngineUtils::LogPclusterData(pcIdGlobal - simulation.pclusters.offset, pidInPcluster, step, data_logging_interval, pos, fe.potE, fe.force, speed,
		simulation.pclusters.count * PersistentCluster::maxParticles, trajBuffer + simulation.logOffset, potEBuffer + simulation.logOffset,
		velocityBuffer + simulation.logOffset, forceBuffer + simulation.logOffset);

	superClusters[scIdGlobal].posX[threadIdx.x] = pos.x;
	superClusters[scIdGlobal].posY[threadIdx.x] = pos.y;
	superClusters[scIdGlobal].posZ[threadIdx.x] = pos.z;
	pclusters[pcIdGlobal].pqd[pidInPcluster].position = pos;
}


