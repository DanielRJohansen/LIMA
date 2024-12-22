//#pragma once - this file must NOT be included multiple times
#pragma once

#include <cuda_runtime.h>


#include "LimaTypes.cuh"
#include "Constants.h"
#include "Simulation.cuh"
#include "LimaPositionSystem.cuh"
#include "KernelConstants.cuh"
#include "BoundaryCondition.cuh"
#include "PhysicsUtilsDevice.cuh"




namespace Electrostatics {





	__device__ inline int convertRelative3dIndexToAbsolute1d(const NodeIndex& nodeIndex) {
		return (nodeIndex.x + 1) + (nodeIndex.y+1) * 3 + (nodeIndex.z+1) * 3*3;
	}

	__device__ inline NodeIndex ConvertAbsolute1dToRelative3d(int index) {
		NodeIndex nodeIndex;

		nodeIndex.z = index / 9 - 1;
		index = index % 9;
		nodeIndex.y = index / 3 - 1;
		nodeIndex.x = index % 3 - 1;

		return nodeIndex;
	}

	__device__ inline NodeIndex GetParticlesBlockRelativeToCompoundOrigo(const Float3& relpos, bool particleActive) {
		if (particleActive) {
			NodeIndex relativeLocalIndex;
			relativeLocalIndex = LIMAPOSITIONSYSTEM::PositionToNodeIndex(relpos);
			relativeLocalIndex.x = std::max(relativeLocalIndex.x, -1);
			relativeLocalIndex.x = std::min(relativeLocalIndex.x, 1);
			relativeLocalIndex.y = std::max(relativeLocalIndex.y, -1);
			relativeLocalIndex.y = std::min(relativeLocalIndex.y, 1);
			relativeLocalIndex.z = std::max(relativeLocalIndex.z, -1);
			relativeLocalIndex.z = std::min(relativeLocalIndex.z, 1);
			return relativeLocalIndex;
		}
		return NodeIndex{};
	}

	// Utilitybuffer min size = sizeof(int) * (27 * 2 + MAX_COMPOUND_PARTICLES)
	//__device__ static void DistributeChargesToChargegrid(const NodeIndex& compoundOrigo, const Float3& relposNM, float charge, ChargeNode* chargeGrid, int nParticles, char* utilityBuffer_sharedMem) {
	__global__ void DistributeCompoundchargesToGridKernel(SimulationDevice * sim) 
	{
		__shared__ int numParticlesInNodeLocal[27];
		__shared__ int numParticlesInNodeGlobal[27];

		NodeIndex compoundOrigo = sim->boxState->compoundOrigos[blockIdx.x];
		const Float3 relPos = sim->boxState->compoundsRelposNm[blockIdx.x * MAX_COMPOUND_PARTICLES];
		const int nParticles = sim->boxConfig.compounds[blockIdx.x].n_particles;
		const float myCharge = sim->boxConfig.compoundsAtomCharges[blockIdx.x * MAX_COMPOUND_PARTICLES + threadIdx.x];

		// todo: threads with a charge of 0 should NOT push that charge (duh)

		if (threadIdx.x < 27) {
			numParticlesInNodeLocal[threadIdx.x] = 0;
			numParticlesInNodeGlobal[threadIdx.x] = 0;
		}
		__syncthreads();

		// First each particle figure out which node it belongs to relative to it's compounds origo. Then it decide the offset to put values, wrt other particles in this compound
		int myOffsetInLocalNode = 0;																				// Between 0 and maxParticlesInNode
		const NodeIndex relativeNode = GetParticlesBlockRelativeToCompoundOrigo(relPos, threadIdx.x < nParticles);	// Between {-1,-1,-1} and {1,1,1}
		const int localNodeId = convertRelative3dIndexToAbsolute1d(relativeNode);									// Between 0 and 26	
		for (int i = 0; i < nParticles; i++) {
			if (threadIdx.x == i)
				myOffsetInLocalNode = numParticlesInNodeLocal[localNodeId]++;
			__syncthreads();
		}

		// Fetch the current particles in the global node counts, and figure out where to push our counts
		if (threadIdx.x < 27) {
			NodeIndex absIndex3d = compoundOrigo + ConvertAbsolute1dToRelative3d(threadIdx.x);
			PeriodicBoundaryCondition::applyBC(absIndex3d);

			// No need to fetch this data, if no particles goes there...
			if (numParticlesInNodeLocal[threadIdx.x] > 0) {
				numParticlesInNodeGlobal[threadIdx.x] = BoxGrid::GetNodePtr(sim->chargeGrid, absIndex3d)->MakeReservation(blockIdx.x, numParticlesInNodeLocal[threadIdx.x]);
			} 
		}
		__syncthreads();

		// Finally compute the correct index to insert our data, and push that data
		if (threadIdx.x < nParticles) {
			NodeIndex absoluteTargetIndex = compoundOrigo + relativeNode;
			PeriodicBoundaryCondition::applyBC(absoluteTargetIndex);

			const int offset = numParticlesInNodeGlobal[convertRelative3dIndexToAbsolute1d(relativeNode)] + myOffsetInLocalNode;
	/*		if (offset < 0 || offset >= ChargeNode::maxParticlesInNode) {
				printf("Error: offset out of bounds: %d. LO %d GO %d\n", offset, localOffset, numParticlesInNodeGlobal[convertRelative3dIndexToAbsolute1d(relativeLocalIndex)]);
				LIMAPOSITIONSYSTEM::PositionToNodeIndex(relPos).print('n');
				(relPos).print('p');
			}*/
			
			//const Float3 positionRelativeToNodeNM = relPos - relativeLocalIndex.toFloat3();
			
			BoxGrid::GetNodePtr(sim->chargeGrid, absoluteTargetIndex)->charges[offset] = myCharge;

			auto res = BoxGrid::GetNodePtr(sim->chargeGrid, absoluteTargetIndex)->compoundReservations[0];
			//BoxGrid::GetNodePtr(sim->chargeGrid, absoluteTargetIndex)->positions[offset] = positionRelativeToNodeNM;
			//BoxGrid::GetNodePtr(sim->chargeGrid, absoluteTargetIndex)->compoundIds[offset] = blockIdx.x;
			//BoxGrid::GetNodePtr(sim->chargeGrid, absoluteTargetIndex)->particleIds[offset] = threadIdx.x;
		}
	}


	//
	__device__ void BitonicSort(ChargeNode::CompoundReservation* reservations) {
		static const int arraySize = ChargeNode::maxReservations;
		for (int k = 2; k <= arraySize; k *= 2) {
			for (int j = k / 2; j > 0; j /= 2) {
				int ixj = threadIdx.x ^ j;
				if (ixj > threadIdx.x) {
					if ((threadIdx.x & k) == 0) {
						if (reservations[threadIdx.x].compoundId > reservations[ixj].compoundId) {
							// Swap
							ChargeNode::CompoundReservation temp = reservations[threadIdx.x];
							reservations[threadIdx.x] = reservations[ixj];
							reservations[ixj] = temp;
						}
					}
					else {
						if (reservations[threadIdx.x].compoundId < reservations[ixj].compoundId) {
							// Swap
							ChargeNode::CompoundReservation temp = reservations[threadIdx.x];
							reservations[threadIdx.x] = reservations[ixj];
							reservations[ixj] = temp;
						}
					}
				}
				__syncthreads();
			}
		}
	}

	static const int SumChargesInGridnode_NTHREADS = ChargeNode::maxReservations;
	__global__ static void SumChargesInGridnode(SimulationDevice* simDev) {
		//__shared__ float chargesBuffer[ChargeNode::maxParticlesInNode];	// OPTIM: each thread can do a parallel partial sum, and then we only need to sum the 32 partial sums in the end

		//ChargeNode* myNode_GlobalMem = BoxGrid::GetNodePtr(simDev->chargeGrid, blockIdx.x);
		//for (int i = threadIdx.x; i < ChargeNode::maxParticlesInNode; i+=blockDim.x) {
		//	chargesBuffer[i] = i < myNode_GlobalMem->nParticles
		//		? myNode_GlobalMem->charges[i]
		//		: 0.f;
		//}
		//__syncthreads();

		//LAL::distributedSummation(chargesBuffer, ChargeNode::maxParticlesInNode);

		//__syncthreads();
		//if (threadIdx.x == 0) {
		//	*BoxGrid::GetNodePtr(simDev->chargeGridChargeSums, blockIdx.x) = chargesBuffer[0];
		//}

		
		__shared__ ChargeNode::CompoundReservation reservations[ChargeNode::maxReservations];
		reservations[threadIdx.x] = ChargeNode::CompoundReservation{};
		//__shared__ float chargesBuffer[ChargeNode::maxParticlesInNode];	// OPTIM: each thread can do a parallel partial sum, and then we only need to sum the 32 partial sums in the end
		//__shared__ float chargesSums[ChargeNode::maxReservations];
		__shared__ float chargeSum;

		ChargeNode* myNode_GlobalMem = BoxGrid::GetNodePtr(simDev->chargeGrid, blockIdx.x);
		const uint32_t nReservations = myNode_GlobalMem->reservationKey >> 16;
		if (threadIdx.x < nReservations)
			reservations[threadIdx.x] = myNode_GlobalMem->compoundReservations[threadIdx.x];
		//chargesSums[threadIdx.x] = 0.f;
		chargeSum = 0.f;
		__syncthreads();

		BitonicSort(reservations);
		__syncthreads();



		float chargeSumInterrim = 0.f;
		if (threadIdx.x < nReservations) {
			for (int i = 0; i < reservations[threadIdx.x].nParticlesInThisNode; i++) {
				chargeSumInterrim += myNode_GlobalMem->charges[reservations[threadIdx.x].firstParticlesOffsetInThisNode + i];
			}
		}

		for (int i = 0; i < nReservations; i++) {
			if (threadIdx.x == i)
				chargeSum += chargeSumInterrim;
			__syncthreads;
		}

		if (threadIdx.x == 0) {
			const NodeIndex mo = BoxGrid::Get3dIndex(blockIdx.x);

			//myNode_GlobalMem->charges[0] = chargeSum;
			myNode_GlobalMem->reservationKey = 0;
			*BoxGrid::GetNodePtr(simDev->chargeGridChargeSums, blockIdx.x) = chargeSum;
		}
	}

	const int CalcLongrangeElectrostaticForces_nThreads = 64;
	__global__ static void CalcLongrangeElectrostaticForces(SimulationDevice* simDev) {
		__shared__ Float3 forceInterims[CalcLongrangeElectrostaticForces_nThreads];	// Each thread accumulates forces from the nodes it has seen
		__shared__ float potEInterims[CalcLongrangeElectrostaticForces_nThreads];
		forceInterims[threadIdx.x] = Float3{};
		potEInterims[threadIdx.x] = 0.f;

		const NodeIndex myNodeindex = BoxGrid::Get3dIndex(blockIdx.x);
		//if (BoxGrid::GetNodePtr(simDev->chargeGridOutputForceAndPot, myNodeindex) == nullptr)
		//	printf("nullptr 1");
		//if (BoxGrid::GetNodePtr(simDev->chargeGrid, myNodeindex)->nParticles == 0)
		//	return;


		// For even sized grids, there is 1 node that is equally far away on both sides, thus that nodes forces cancel out. This logic avoids such nodes
		const int nNodesPerDimRelativeToUs = LAL::IsEven(boxSize_device.blocksPerDim) ? boxSize_device.blocksPerDim - 1 : boxSize_device.blocksPerDim;
		const int nNodesInBoxGridEffective = BoxGrid::BlocksTotal(nNodesPerDimRelativeToUs);
		const NodeIndex nodeOffset{ -nNodesPerDimRelativeToUs / 2,-nNodesPerDimRelativeToUs / 2,-nNodesPerDimRelativeToUs / 2 };

		for (int i = threadIdx.x; i < nNodesInBoxGridEffective; i += blockDim.x) {
			const NodeIndex queryNodeindexRelative = BoxGrid::Get3dIndexWithNNodes(i, nNodesPerDimRelativeToUs) + nodeOffset;

			if (queryNodeindexRelative.MaxAbsElement() < 2) // These are shortrange 
				continue;

			NodeIndex queryNodeindexAbsolute = myNodeindex + queryNodeindexRelative;
			PeriodicBoundaryCondition::applyBC(queryNodeindexAbsolute);

			if (BoxGrid::GetNodePtr(simDev->chargeGridOutputForceAndPot, queryNodeindexAbsolute) == nullptr)
				printf("nullptr 2");

			const Float3 diff = -queryNodeindexRelative.toFloat3() * static_cast<float>(BoxGrid::blocksizeNM); // [nm]
			const float queryCharge = *BoxGrid::GetNodePtr(simDev->chargeGridChargeSums, queryNodeindexAbsolute);

			forceInterims[threadIdx.x] += PhysicsUtilsDevice::CalcCoulumbForce_optim(1.f * queryCharge, diff) * PhysicsUtilsDevice::modifiedCoulombConstant_Force;
			potEInterims[threadIdx.x] += PhysicsUtilsDevice::CalcCoulumbPotential_optim(1.f, queryCharge, diff) * 0.5f * PhysicsUtilsDevice::modifiedCoulombConstant_Potential;	// 0.5 because the other node is also calcing this
		}

		__syncthreads();
		LAL::distributedSummation(forceInterims, CalcLongrangeElectrostaticForces_nThreads);
		LAL::distributedSummation(potEInterims, CalcLongrangeElectrostaticForces_nThreads);
		__syncthreads();

		if (threadIdx.x == 0) {
			// TODO: this potential should be split between all particles in the chargenode, RIIIIIIIIIIIIIIIIGHT?? Josiah??
			//printf("Force out: %f pot out %f\n", forceInterims[0].len(), potEInterims[0]);

			
			BoxGrid::GetNodePtr(simDev->chargeGridOutputForceAndPot, myNodeindex)->forcePart = forceInterims[0];
			BoxGrid::GetNodePtr(simDev->chargeGridOutputForceAndPot, myNodeindex)->potentialPart = potEInterims[0];			
		}
	}


	// Returns timing in [ys]
	__host__ static int HandleElectrostatics(SimulationDevice* sim_dev, BoxParams boxparamsHost, cudaStream_t& cudaStream) 
	{
		// TODO Needs a cudastream as argument
		//LIMA_UTILS::genericErrorCheckNoSync("Error Before Electrostatics SumChargesInGridnode");
		const auto t0 = std::chrono::high_resolution_clock::now();

		//static_assert(ChargeNode::maxParticlesInNode %32 == 0, "Chargenode charges size isn't optimal for this kernel");
		SumChargesInGridnode<<<BoxGrid::BlocksTotal(boxparamsHost.boxSize), SumChargesInGridnode_NTHREADS, 0, cudaStream >>>(sim_dev);
		//LIMA_UTILS::genericErrorCheck("Error after Electrostatics SumChargesInGridnode");


		CalcLongrangeElectrostaticForces<<<BoxGrid::BlocksTotal(boxparamsHost.boxSize), CalcLongrangeElectrostaticForces_nThreads, 0, cudaStream>>>(sim_dev);
		//LIMA_UTILS::genericErrorCheck("Error after Electrostatics CalcLongrangeElectrostaticForces");


		const auto t1 = std::chrono::high_resolution_clock::now();
		return static_cast<int>(std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count());
	}


} // namespace Electrostatics