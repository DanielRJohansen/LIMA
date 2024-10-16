#pragma once


#include "LimaTypes.cuh"
#include "Constants.h"
#include "Simulation.cuh"
#include "EngineUtilsWarnings.cuh"
#include "BoundaryConditionPublic.h"
#include "LimaPositionSystem.cuh"
#include "LimaTypes.cuh"
#include "Constants.h"
#include "Bodies.cuh"
#include "BoxGrid.cuh"
#include "PhysicsUtils.cuh"
#include "SimulationDevice.cuh"
#include "KernelWarnings.cuh"
#include "KernelConstants.cuh"

#include <cooperative_groups.h>
#include <cooperative_groups/memcpy_async.h>
#include <curand_kernel.h>  // For generating random numbers

namespace EngineUtils {

	template <typename BoundaryCondition>
	__device__ int static getNewBlockId(const NodeIndex& transfer_direction, const NodeIndex& origo) {
		NodeIndex new_nodeindex = transfer_direction + origo;
		BoundaryCondition::applyBC(new_nodeindex);
		return BoxGrid::Get1dIndex(new_nodeindex, boxSize_device.boxSizeNM_i);
	}

	// returns pos_tadd1
	__device__ static Coord integratePositionVVS(const Coord& pos, const Float3& vel, const Float3& force, const float mass, const float dt) {
#ifndef ENABLE_INTEGRATEPOSITION
		return pos;
#endif
		const Coord pos_tadd1 = pos + Coord{ (vel * dt + force * (0.5 / mass * dt * dt)).round() };				// precise version
		return pos_tadd1;
	}
	__device__ static Float3 integrateVelocityVVS(const Float3& vel_tsub1, const Float3& force_tsub1, const Float3& force, const float dt, const float mass) {
		const Float3 vel = vel_tsub1 + (force + force_tsub1) * (dt * 0.5f / mass);
		return vel;
	}


	// ChatGPT magic. generates a float with elements between -1 and 1
	__device__ inline Float3 GenerateRandomForce() {
		unsigned int seed = threadIdx.x + blockIdx.x*blockDim.x;

		// Simple LCG (Linear Congruential Generator) for pseudo-random numbers
		seed = (1664525 * seed + 1013904223);
		float randX = ((seed & 0xFFFF) / 32768.0f) - 1.0f;

		seed = (1664525 * seed + 1013904223);
		float randY = ((seed & 0xFFFF) / 32768.0f) - 1.0f;

		seed = (1664525 * seed + 1013904223);
		float randZ = ((seed & 0xFFFF) / 32768.0f) - 1.0f;

		return Float3{randX, randY, randZ};
	}

	// Tanh activation functions that scales forces during EM
	__device__ static Float3 ForceActivationFunction(const Float3 force, float progress /*How far along the EM are we, 0..1*/) {

		// Handled inf forces by returning a pseudorandom z force based on global thread index
		if (isinf(force.lenSquared())) {
			return GenerateRandomForce();
		}


		// 1000 [kJ/mol/nm] is a good target for EM. For EM we will scale the forces below this value * 200
		const float scaleAbove = 1000.f + 30000.f * (1.f-progress);
		const float alpha = scaleAbove * LIMA / NANO * KILO; // [1/l N/mol]

		// Apply tanh to the magnitude

		const float scaledMagnitude = alpha * tanh(force.len()/alpha);
		// Scale the original force vector by the ratio of the new magnitude to the original magnitude
		Float3 scaledForce = force * (scaledMagnitude / (force.len() + 1e-8)); // Avoid division by zero		
		return scaledForce;
	}

	/*__device__ static Float3 SlowHighEnergyParticle(const Float3& velocityNow, const float dt, const float mass, float thermostatScalar) {
		const float kineticEnergy = PhysicsUtils::calcKineticEnergy(velocityNow.len(), mass);
		const float temperature = PhysicsUtils::kineticEnergyToTemperature(kineticEnergy);

		const float temperatureThreshold = 400.f; // [K]

		if (temperature > temperatureThreshold) {
			return velocityNow * 0.9995f;
		}
		else {
			return velocityNow * thermostatScalar;
		}
	}*/





	__device__ inline void LogCompoundData(const CompoundCompact& compound, int totalParticlesUpperbound, CompoundCoords& compound_coords, 
		float* potE_sum, const Float3& force, Float3& force_LJ_sol, const SimParams& simparams, SimSignals& simsignals, 
		float* poteBuffer, Float3* trajBuffer, float* velBuffer, Float3* forceBuffer, const float speed, int64_t step)
	{
		if (threadIdx.x >= compound.n_particles) { return; }

		if (step % simparams.data_logging_interval != 0) { return; }

		const int index = DatabuffersDeviceController::GetLogIndexOfParticle(threadIdx.x, blockIdx.x, step, simparams.data_logging_interval, totalParticlesUpperbound);
		trajBuffer[index] = LIMAPOSITIONSYSTEM::GetAbsolutePositionNM(compound_coords.origo, compound_coords.rel_positions[threadIdx.x]); 
		poteBuffer[index] = *potE_sum;
		velBuffer[index] = speed;
		forceBuffer[index] = force;

		EngineUtilsWarnings::logcompoundVerifyVelocity(compound, simparams, simsignals, compound_coords, force, speed);
	}

	__device__ inline void LogSolventData(const BoxParams& boxparams, const float& potE, const SolventBlock& solventblock, bool solvent_active, 
		const Float3& force, const Float3& velocity, uint32_t step, float* poteBuffer, Float3* trajBuffer, float* velBuffer, int loggingInterval)
	{
		if (step % loggingInterval != 0) { return; }


		if (solvent_active) {
			const int index = DatabuffersDeviceController::GetLogIndexOfParticle(solventblock.ids[threadIdx.x], boxparams.n_compounds, step, 
				loggingInterval, boxparams.total_particles_upperbound);

			trajBuffer[index] = LIMAPOSITIONSYSTEM::GetAbsolutePositionNM(solventblock.origo, solventblock.rel_pos[threadIdx.x]);
			poteBuffer[index] = potE;
			velBuffer[index] = velocity.len();
		}
	}

	__device__ constexpr bool isOutsideCutoff(const float dist_sq_reciprocal) {
		if constexpr (HARD_CUTOFF) {
			return dist_sq_reciprocal < cutoffLmSquaredReciprocal_device;	//  1. / (CUTOFF_LM * CUTOFF_LM);
		}
		return false;
	}

	//template <bool useDefault>
	//__device__ constexpr bool isOutsideCutoff(const float dist_sq_reciprocal);

	//template <>
	//__device__ constexpr bool isOutsideCutoff<true>(const float dist_sq_reciprocal) {
	//	if constexpr (HARD_CUTOFF) {
	//		const float threshold = cutoffLmSquaredReciprocal_device;
	//		return dist_sq_reciprocal < 1.2f;	//  1. / (CUTOFF_LM * CUTOFF_LM);
	//	}
	//	return false;
	//}
	//template <>
	//__device__ constexpr bool isOutsideCutoff<false>(const float dist_sq_reciprocal) {
	//	if constexpr (HARD_CUTOFF) {			
	//		return dist_sq_reciprocal < cutoffLmSquaredReciprocal_device; //  1. / (CUTOFF_LM * CUTOFF_LM);
	//	}
	//	return false;
	//}









	template <typename BoundaryCondition>
	__device__ inline void getCompoundHyperpositionsAsFloat3(const NodeIndex& origo_self, const CompoundCoords* const querycompound,
		void* output_buffer, Float3& utility_float3, const int n_particles)
	{
		if (threadIdx.x == 0) {
			const NodeIndex querycompound_hyperorigo = BoundaryCondition::applyHyperpos_Return(origo_self, querycompound->origo);
			KernelHelpersWarnings::assertHyperorigoIsValid(querycompound_hyperorigo, origo_self);

			// calc Relative LimaPosition Shift from the origo-shift
			utility_float3 = LIMAPOSITIONSYSTEM_HACK::getRelShiftFromOrigoShift(querycompound_hyperorigo, origo_self).toFloat3();
		}
		//	__syncthreads();

		auto block = cooperative_groups::this_thread_block();
		cooperative_groups::memcpy_async(block, (Coord*)output_buffer, querycompound->rel_positions, sizeof(Coord) * n_particles);
		cooperative_groups::wait(block);
		__syncthreads();

		// Eventually i could make it so i only copy the active particles in the compound
		if (threadIdx.x < n_particles) {
			const Coord queryparticle_coord = ((Coord*)output_buffer)[threadIdx.x];
			((Float3*)output_buffer)[threadIdx.x] = queryparticle_coord.toFloat3() + utility_float3;
		}
		__syncthreads();
	}

	__device__ inline void getCompoundHyperpositionsAsFloat3Async(const CompoundCoords* const querycompound,
		void* output_buffer, const int n_particles, const Float3& relshift)
	{
		auto block = cooperative_groups::this_thread_block();
		cooperative_groups::memcpy_async(block, (Coord*)output_buffer, querycompound->rel_positions, sizeof(Coord) * n_particles);
		cooperative_groups::wait(block);

		// Eventually i could make it so i only copy the active particles in the compound
		if (threadIdx.x < n_particles) {
			static_assert(sizeof(Float3) == sizeof(Coord), "Float3 and Coord must have same size");
			/*const Coord queryparticle_coord = ((Coord*)output_buffer)[threadIdx.x];
			((Float3*)output_buffer)[threadIdx.x] = queryparticle_coord.toFloat3() + relshift;*/
		}
	}

	//__device__ inline void LoadNextTaskAsync(int taskIndex, int indexInBatch, const NeighborList& neighborlist, 
	//	int neighbor_n_particles[32], uint8_t* neighborAtomstypesNext, BoxDevice* box, thread_block block) {
	//	//static_assert(sizeof(Coord) == sizeof(Float3));
	//	const int nextNeighborIndex = neighborlist.neighborcompound_ids[taskIndex + 1];
	//	const int nextNeighborNParticles = neighbor_n_particles[indexInBatch + 1];
	//	neighborAtomstypesNext[threadIdx.x] = box->compounds[nextNeighborIndex].atom_types[threadIdx.x];// TODO figure out how to make this async. Probably by having all atomstypes in a single global buffer...
	//	//neighborParticleschargesNext[threadIdx.x] = box->compounds[nextNeighborIndex].atom_charges[threadIdx.x];

	//	cooperative_groups::memcpy_async(block, (Coord*)neighborPositionsNext, coords_ptrs[indexInBatch + 1]->rel_positions, sizeof(Coord) * nextNeighborNParticles);
	//	cooperative_groups::memcpy_async(block, neighborParticleschargesNext, box->compounds[nextNeighborIndex].atom_charges, sizeof(half) * nextNeighborNParticles);
	//	if (i + 1 >= compound.n_bonded_compounds)
	//		cooperative_groups::memcpy_async(block, bplutNext, (BondedParticlesLUT*)compoundPairLUTs[indexInBatch + 1], sizeof(BondedParticlesLUT));
	//	//cooperative_groups::memcpy_async(block, neighborAtomstypesNext, (uint8_t*)box->compounds[nextNeighborIndex].atom_types, sizeof(uint8_t)* MAX_COMPOUND_PARTICLES);				

	//}

	

};


namespace LIMAKERNELDEBUG {

	//__device__ void static compoundIntegration(const Coord& relpos_prev, const Coord& relpos_next, const Float3& force, bool& critical_error_encountered) {
	//	const auto dif = (relpos_next - relpos_prev);
	//	const int32_t max_diff = BOXGRID_NODE_LEN_i / 20;
	//	if (std::abs(dif.x) > max_diff || std::abs(dif.y) > max_diff || std::abs(dif.z) > max_diff || force.len() > 3.f) {
	//		//printf("\nParticle %d in compound %d is moving too fast\n", threadIdx.x, blockIdx.x);
	//		printf("\nParticle is moving too fast\n");
	//		dif.printS('D');
	//		force.print('F');
	//		relpos_next.printS('R');
	//		critical_error_encountered = true;
	//	}
	//}

	//__device__ void static solventIntegration(const Coord& relpos_prev, const Coord& relpos_next, const Float3& force, bool& critical_error_encountered, int id) {
	//	const auto dif = (relpos_next - relpos_prev);
	//	const int32_t max_diff = BOXGRID_NODE_LEN_i / 20;
	//	if (std::abs(dif.x) > max_diff || std::abs(dif.y) > max_diff || std::abs(dif.z) > max_diff) {
	//		printf("\nSolvent %d moving too fast\n", id);
	//		dif.printS('D');
	//		force.print('F');
	//		critical_error_encountered = true;
	//	}
	//}

	// 
	//__device__ bool static forceTooLarge(const Float3& force, const float threshold=0.5) { // 

	//}

};

namespace DEBUGUTILS {

	/// <summary>
	/// Puts the nearest solvent of each solvent in the out vector.
	/// </summary>
	//void findAllNearestSolventSolvent(SolventBlocksCircularQueue* queue, size_t n_solvents, std::vector<float>& out);
}

