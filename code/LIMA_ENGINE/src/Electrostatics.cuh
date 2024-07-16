//#pragma once - this file must NOT be included multiple times
#pragma once

#include <cuda_runtime.h>


#include "LimaTypes.cuh"
#include "Constants.h"
#include "Simulation.cuh"
#include "LimaPositionSystem.cuh"
#include "KernelConstants.cuh"
#include "BoundaryCondition.cuh"










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

	__device__ static void DistributeChargesToChargegrid(const NodeIndex& compoundOrigo, const Float3& relposLM, float charge, ChargeNode* chargeGrid, int nParticles, char* utilityBuffer_sharedMem) {
		

		int* numParticlesInNodeLocal = (int*)utilityBuffer_sharedMem;

		// First clean the memory
		static_assert(MAX_COMPOUND_PARTICLES >= 27 * 2, "Not enough threads to reset buffer");
		if (threadIdx.x < 27 * 2) {
			numParticlesInNodeLocal[threadIdx.x] = 0;
		}
		__syncthreads();


		int localOffset = 0;
		NodeIndex relativeLocalIndex;
		if (threadIdx.x < nParticles) {
			relativeLocalIndex = LIMAPOSITIONSYSTEM::PositionToNodeIndex(relposLM);
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

		int* numParticlesInNodeGlobal = &((int*)utilityBuffer_sharedMem)[27];
		if (threadIdx.x < 27) {
			NodeIndex absIndex3d = compoundOrigo + ConvertAbsolute1dToRelative3d(threadIdx.x);
			PeriodicBoundaryCondition::applyBC(absIndex3d);	

			int* nParticlesInNodePtr = &BoxGrid::GetNodePtr(chargeGrid, absIndex3d)->nParticles;
			numParticlesInNodeGlobal[threadIdx.x] = atomicAdd(nParticlesInNodePtr, numParticlesInNodeLocal[threadIdx.x]);
		}
		__syncthreads();

		if (threadIdx.x < nParticles) {
			NodeIndex absoluteTargetIndex = compoundOrigo + relativeLocalIndex;
			PeriodicBoundaryCondition::applyBC(absoluteTargetIndex);

			const int offset = numParticlesInNodeGlobal[convertRelative3dIndexToAbsolute1d(relativeLocalIndex)] + localOffset;
			if (offset < 0 || offset >= ChargeNode::maxParticlesInNode) {
				printf("Error: offset out of bounds: %d. LO %d GO %d\n", offset, localOffset, numParticlesInNodeGlobal[convertRelative3dIndexToAbsolute1d(relativeLocalIndex)]);
				LIMAPOSITIONSYSTEM::PositionToNodeIndex(relposLM).print('n');
				relposLM.print('r');
				(relposLM * LIMA_TO_NANO).print('p');
			}
			
			const Float3 positionRelativeToNodeNM = relposLM * LIMA_TO_NANO - relativeLocalIndex.toFloat3();
			BoxGrid::GetNodePtr(chargeGrid, absoluteTargetIndex)->positions[offset] = positionRelativeToNodeNM;
			BoxGrid::GetNodePtr(chargeGrid, absoluteTargetIndex)->charges[offset] = charge;
			BoxGrid::GetNodePtr(chargeGrid, absoluteTargetIndex)->compoundIds[offset] = blockIdx.x;
			BoxGrid::GetNodePtr(chargeGrid, absoluteTargetIndex)->particleIds[offset] = threadIdx.x;
		}
	}






	__global__ static void HandleShortrangeElectrostatics(SimulationDevice* simDev) {
		__shared__ Float3 positionsBuffer[ChargeNode::maxParticlesInNode];
		__shared__ float chargesBuffer[ChargeNode::maxParticlesInNode];

		const int index1D = blockIdx.x;
		const NodeIndex index3D = BoxGrid::Get3dIndex(blockIdx.x);
		ChargeNode* myNode_GlobalMem = BoxGrid::GetNodePtr(simDev->chargeGrid, index1D);
		const bool threadActive = threadIdx.x < myNode_GlobalMem->nParticles;

		const Float3 myPos = threadActive
			? myNode_GlobalMem->positions[threadIdx.x]
			: Float3(42.f);

		const float myCharge = threadActive
			? myNode_GlobalMem->charges[threadIdx.x]
			: 0.f;

		Float3 force = 0.f; // Accumulated as [GN/mol], converted to GigaN before writing to output
		float potE = 0.f; // Accumulated as [J/mol], converted to GigaJ before writing to output


		for (int xOff = -1; xOff <= 1; xOff++) {
			for (int yOff = -1; yOff <= 1; yOff++) {
				for (int zOff = -1; zOff <= 1; zOff++) {

					const bool skipOwnIndex = xOff == 0 && yOff == 0 && zOff == 0;
					NodeIndex queryNodeindex = index3D + NodeIndex(xOff, yOff, zOff);
					PeriodicBoundaryCondition::applyBC(queryNodeindex);

					if (threadIdx.x < BoxGrid::GetNodePtr(simDev->chargeGrid, queryNodeindex)->nParticles) {
						positionsBuffer[threadIdx.x] = LIMAPOSITIONSYSTEM::GetAbsolutePositionNM(NodeIndex{xOff, yOff, zOff}, BoxGrid::GetNodePtr(simDev->chargeGrid, queryNodeindex)->positions[threadIdx.x]);
						
						chargesBuffer[threadIdx.x] = BoxGrid::GetNodePtr(simDev->chargeGrid, queryNodeindex)->charges[threadIdx.x];
					}
					__syncthreads();

					if (threadActive) {
						for (int i = 0; i < BoxGrid::GetNodePtr(simDev->chargeGrid, queryNodeindex)->nParticles; i++) {

							if (skipOwnIndex && i == threadIdx.x)
								continue;

							const Float3 otherPos = positionsBuffer[i];
							const float otherCharge = chargesBuffer[i];

							const Float3 diff = myPos - otherPos;	
								
							force += PhysicsUtils::CalcCoulumbForce(myCharge, otherCharge, diff);
							potE += PhysicsUtils::CalcCoulumbPotential(myCharge, otherCharge, diff.len()) * 0.5f;	// 0.5 because the other particle is also calcing this
						}
					}
				}
			}
		}




		if (threadActive) {
			const int cid = myNode_GlobalMem->compoundIds[threadIdx.x];
			const int pid = myNode_GlobalMem->particleIds[threadIdx.x];
			//printf("Resulting force %f potE %f \n", force.len(), potE);

			simDev->box->compounds[cid].potE_interim[pid] += potE;
			simDev->box->compounds[cid].forces_interim[pid] += force; 
		}

		__syncthreads();

		// Finally reset nParticles before next step
		if (threadIdx.x == 0)
			BoxGrid::GetNodePtr(simDev->chargeGrid, index3D)->nParticles = 0;
	}

	__host__ static void HandleElectrostatics(SimulationDevice* sim_dev, BoxParams boxparamsHost) 
	{
		// First handle short range
		HandleShortrangeElectrostatics<<<BoxGrid::BlocksTotal(boxparamsHost.boxSize), ChargeNode::maxParticlesInNode >> > (sim_dev);
		//HandleShortrangeElectrostatics << <32, 7 >> > (sim_dev);
	}


} // namespace Electrostatics