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

	// Utilitybuffer min size = sizeof(int) * (27 * 2 + MAX_COMPOUND_PARTICLES)
	__device__ static void DistributeChargesToChargegrid(const NodeIndex& compoundOrigo, const Float3& relposNM, float charge, ChargeNode* chargeGrid, int nParticles, char* utilityBuffer_sharedMem) {
		int* numParticlesInNodeLocal = (int*)utilityBuffer_sharedMem; // First 27 ints

		// First clean the memory
		//static_assert(MAX_COMPOUND_PARTICLES >= 27 * 2, "Not enough threads to reset buffer");
		//if (threadIdx.x < 27 * 2) {
		for (int i = threadIdx.x; i < 27 * 2; i += blockDim.x) {
			numParticlesInNodeLocal[threadIdx.x] = 0;
		}
		__syncthreads();

		// First each particle figure out which node it belongs to RELATIVE TO ITS COMPOUNDS ORIGO
		int localOffset = 0;
		NodeIndex relativeLocalIndex;
		if (threadIdx.x < nParticles) {
			relativeLocalIndex = LIMAPOSITIONSYSTEM::PositionToNodeIndex(relposNM);
			relativeLocalIndex.x = std::max(relativeLocalIndex.x, -1);
			relativeLocalIndex.x = std::min(relativeLocalIndex.x, 1);
			relativeLocalIndex.y = std::max(relativeLocalIndex.y, -1);
			relativeLocalIndex.y = std::min(relativeLocalIndex.y, 1);
			relativeLocalIndex.z = std::max(relativeLocalIndex.z, -1);
			relativeLocalIndex.z = std::min(relativeLocalIndex.z, 1);


			const int localIndex = convertRelative3dIndexToAbsolute1d(relativeLocalIndex);
			localOffset = atomicAdd(&numParticlesInNodeLocal[localIndex], 1);
		}
		__syncthreads();

		// Fetch the current particles in the global node counts, and figure out where to push our counts
		int* numParticlesInNodeGlobal = &((int*)utilityBuffer_sharedMem)[27];
		if (threadIdx.x < 27) {
			NodeIndex absIndex3d = compoundOrigo + ConvertAbsolute1dToRelative3d(threadIdx.x);
			PeriodicBoundaryCondition::applyBC(absIndex3d);

			//if (BoxGrid::GetNodePtr(chargeGrid, absIndex3d) == nullptr)
			//	printf("nullptr 6");

			int* nParticlesInNodePtr = &BoxGrid::GetNodePtr(chargeGrid, absIndex3d)->nParticles;

			numParticlesInNodeGlobal[threadIdx.x] = atomicAdd(nParticlesInNodePtr, numParticlesInNodeLocal[threadIdx.x]);
		/*	if (numParticlesInNodeGlobal[threadIdx.x] > ChargeNode::maxParticlesInNode) {
				printf("Error: Too many particles in node from other kernels: %d\n", numParticlesInNodeGlobal[threadIdx.x]);
			}*/
		}
		__syncthreads();

		// Finally compute the correct index to insert our data, and push that data
		if (threadIdx.x < nParticles) {
			NodeIndex absoluteTargetIndex = compoundOrigo + relativeLocalIndex;
			PeriodicBoundaryCondition::applyBC(absoluteTargetIndex);

		/*	if (BoxGrid::GetNodePtr(chargeGrid, absoluteTargetIndex) == nullptr)
				printf("nullptr 5");*/


			const int offset = numParticlesInNodeGlobal[convertRelative3dIndexToAbsolute1d(relativeLocalIndex)] + localOffset;
			if (offset < 0 || offset >= ChargeNode::maxParticlesInNode) {
				printf("Error: offset out of bounds: %d. LO %d GO %d\n", offset, localOffset, numParticlesInNodeGlobal[convertRelative3dIndexToAbsolute1d(relativeLocalIndex)]);
				LIMAPOSITIONSYSTEM::PositionToNodeIndex(relposNM).print('n');
				(relposNM).print('p');
			}
			
			const Float3 positionRelativeToNodeNM = relposNM - relativeLocalIndex.toFloat3();
			BoxGrid::GetNodePtr(chargeGrid, absoluteTargetIndex)->positions[offset] = positionRelativeToNodeNM;
			BoxGrid::GetNodePtr(chargeGrid, absoluteTargetIndex)->charges[offset] = charge;
			BoxGrid::GetNodePtr(chargeGrid, absoluteTargetIndex)->compoundIds[offset] = blockIdx.x;
			BoxGrid::GetNodePtr(chargeGrid, absoluteTargetIndex)->particleIds[offset] = threadIdx.x;
		}
	}






	//__global__ static void HandleShortrangeElectrostatics(SimulationDevice* simDev) {
	//	__shared__ Float3 positionsBuffer[ChargeNode::maxParticlesInNode];
	//	__shared__ float chargesBuffer[ChargeNode::maxParticlesInNode];

	//	const int index1D = blockIdx.x;
	//	const NodeIndex index3D = BoxGrid::Get3dIndex(blockIdx.x);
	//	ChargeNode* myNode_GlobalMem = BoxGrid::GetNodePtr(simDev->chargeGrid, index1D);
	//	const bool threadActive = threadIdx.x < myNode_GlobalMem->nParticles;

	//	const Float3 myPos = threadActive
	//		? myNode_GlobalMem->positions[threadIdx.x]
	//		: Float3(42.f);

	//	const float myCharge = threadActive
	//		? myNode_GlobalMem->charges[threadIdx.x]
	//		: 0.f;

	//	Float3 force{}; // Accumulated as [GN/mol], converted to GigaN before writing to output
	//	float potE = 0.f; // Accumulated as [J/mol], converted to GigaJ before writing to output


	//	for (int xOff = -1; xOff <= 1; xOff++) {
	//		for (int yOff = -1; yOff <= 1; yOff++) {
	//			for (int zOff = -1; zOff <= 1; zOff++) {

	//				const bool skipOwnIndex = xOff == 0 && yOff == 0 && zOff == 0;
	//				NodeIndex queryNodeindex = index3D + NodeIndex(xOff, yOff, zOff);
	//				PeriodicBoundaryCondition::applyBC(queryNodeindex);

	//				if (threadIdx.x < BoxGrid::GetNodePtr(simDev->chargeGrid, queryNodeindex)->nParticles) {
	//					positionsBuffer[threadIdx.x] = LIMAPOSITIONSYSTEM::GetAbsolutePositionNM(NodeIndex{xOff, yOff, zOff}, BoxGrid::GetNodePtr(simDev->chargeGrid, queryNodeindex)->positions[threadIdx.x]);
	//					
	//					chargesBuffer[threadIdx.x] = BoxGrid::GetNodePtr(simDev->chargeGrid, queryNodeindex)->charges[threadIdx.x];
	//				}
	//				__syncthreads();

	//				if (threadActive) {
	//					for (int i = 0; i < BoxGrid::GetNodePtr(simDev->chargeGrid, queryNodeindex)->nParticles; i++) {

	//						if (skipOwnIndex && i == threadIdx.x)
	//							continue;

	//						const Float3 otherPos = positionsBuffer[i];
	//						const float otherCharge = chargesBuffer[i];

	//						const Float3 diff = myPos - otherPos;	
	//							
	//						force += PhysicsUtilsDevice::CalcCoulumbForce(myCharge, otherCharge, diff);
	//						potE += PhysicsUtilsDevice::CalcCoulumbPotential(myCharge, otherCharge, diff) * 0.5f;	// 0.5 because the other particle is also calcing this
	//					}
	//				}
	//			}
	//		}
	//	}




	//	if (threadActive) {
	//		const int cid = myNode_GlobalMem->compoundIds[threadIdx.x];
	//		const int pid = myNode_GlobalMem->particleIds[threadIdx.x];
	//		//printf("Resulting force %f potE %f \n", force.len(), potE);

	//		simDev->boxState->compounds[cid].potE_interim[pid] += potE;
	//		simDev->boxState->compounds[cid].forces_interim[pid] += force;
	//	}

	//	__syncthreads();

	//	// Finally reset nParticles before next step
	//	if (threadIdx.x == 0)
	//		BoxGrid::GetNodePtr(simDev->chargeGrid, index3D)->nParticles = 0;
	//}



	__global__ static void SumChargesInGridnode(SimulationDevice* simDev) {
		__shared__ float chargesBuffer[ChargeNode::maxParticlesInNode];	// OPTIM: each thread can do a parallel partial sum, and then we only need to sum the 32 partial sums in the end

		ChargeNode* myNode_GlobalMem = BoxGrid::GetNodePtr(simDev->chargeGrid, blockIdx.x);
		for (int i = threadIdx.x; i < ChargeNode::maxParticlesInNode; i+=blockDim.x) {
			chargesBuffer[i] = i < myNode_GlobalMem->nParticles
				? myNode_GlobalMem->charges[i]
				: 0.f;
		}
		__syncthreads();

		LAL::distributedSummation(chargesBuffer, ChargeNode::maxParticlesInNode);

		__syncthreads();
		if (threadIdx.x == 0) {
			*BoxGrid::GetNodePtr(simDev->chargeGridChargeSums, blockIdx.x) = chargesBuffer[0];
		}
	}

	const int CalcLongrangeElectrostaticForces_nThreads = 64;
	__global__ static void CalcLongrangeElectrostaticForces(SimulationDevice* simDev) {
		__shared__ Float3 forceInterims[CalcLongrangeElectrostaticForces_nThreads];	// Each thread accumulates forces from the nodes it has seen
		__shared__ float potEInterims[CalcLongrangeElectrostaticForces_nThreads];
		forceInterims[threadIdx.x] = Float3{};
		potEInterims[threadIdx.x] = 0.f;

		const NodeIndex myNodeindex = BoxGrid::Get3dIndex(blockIdx.x);
		if (BoxGrid::GetNodePtr(simDev->chargeGridOutputForceAndPot, myNodeindex) == nullptr)
			printf("nullptr 1");
		if (BoxGrid::GetNodePtr(simDev->chargeGrid, myNodeindex)->nParticles == 0)
			return;


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

			forceInterims[threadIdx.x] += PhysicsUtilsDevice::CalcCoulumbForce_optim(1.f, queryCharge, diff) * PhysicsUtilsDevice::modifiedCoulombConstant_Force;
			potEInterims[threadIdx.x] += PhysicsUtilsDevice::CalcCoulumbPotential_optim(1.f, queryCharge, diff) * 0.5f * PhysicsUtilsDevice::modifiedCoulombConstant_Potential;	// 0.5 because the other node is also calcing this
		}

		__syncthreads();
		LAL::distributedSummation(forceInterims, CalcLongrangeElectrostaticForces_nThreads);
		LAL::distributedSummation(potEInterims, CalcLongrangeElectrostaticForces_nThreads);
		__syncthreads();

		if (threadIdx.x == 0) {
			// TODO: this potential should be split between all particles in the chargenode, RIIIIIIIIIIIIIIIIGHT?? Josiah??
			//printf("Force out: %f pot out %f\n", forceInterims[0].len(), potEInterims[0]);

			const float nParticles = static_cast<float>(BoxGrid::GetNodePtr(simDev->chargeGrid, myNodeindex)->nParticles);
			BoxGrid::GetNodePtr(simDev->chargeGridOutputForceAndPot, myNodeindex)->forcePart = forceInterims[0];
			BoxGrid::GetNodePtr(simDev->chargeGridOutputForceAndPot, myNodeindex)->potentialPart = potEInterims[0];
			BoxGrid::GetNodePtr(simDev->chargeGrid, myNodeindex)->nParticles = 0;
		}
	}


	// Returns timing in [ys]
	__host__ static int HandleElectrostatics(SimulationDevice* sim_dev, BoxParams boxparamsHost) 
	{
		LIMA_UTILS::genericErrorCheck("Error Before Electrostatics SumChargesInGridnode");
		const auto t0 = std::chrono::high_resolution_clock::now();

		// First handle short range
		//HandleShortrangeElectrostatics<<<BoxGrid::BlocksTotal(boxparamsHost.boxSize), ChargeNode::maxParticlesInNode >> > (sim_dev);

		static_assert(ChargeNode::maxParticlesInNode %32 == 0, "Chargenode charges size isn't optimal for this kernel");
		SumChargesInGridnode<<<BoxGrid::BlocksTotal(boxparamsHost.boxSize), 32>>>(sim_dev);
		LIMA_UTILS::genericErrorCheck("Error after Electrostatics SumChargesInGridnode");


		CalcLongrangeElectrostaticForces<<<BoxGrid::BlocksTotal(boxparamsHost.boxSize), CalcLongrangeElectrostaticForces_nThreads>>>(sim_dev);
		LIMA_UTILS::genericErrorCheck("Error after Electrostatics CalcLongrangeElectrostaticForces");


		const auto t1 = std::chrono::high_resolution_clock::now();
		return static_cast<int>(std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count());
	}


} // namespace Electrostatics