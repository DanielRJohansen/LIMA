#pragma once


#include "LimaTypes.cuh"
#include "Constants.h"
#include "Simulation.cuh"
#include "EngineUtilsWarnings.cuh"
#include "LimaPositionSystem.cuh"
#include "LimaTypes.cuh"
#include "Constants.h"
#include "Bodies.cuh"
#include "BoxGrid.cuh"
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
		const Coord pos_tadd1 = pos + Coord{ (vel * dt + force * (0.5f / mass * dt * dt)) };				// precise version
		return pos_tadd1;
	}
	__device__ static Float3 integrateVelocityVVS(const Float3& vel_tsub1, const Float3& force_tsub1, const Float3& force, const float dt, const float mass) {
		const Float3 vel = vel_tsub1 + (force + force_tsub1) * (dt * 0.5f / mass);
		return vel;
	}
	

	__device__ static Coord IntegratePositionADAM(const Coord& pos, Float3 force, AdamState* const adamState, int step) {

		// TODO: Maybe figure out the highest force particle in the system, and scale each particles lr based on that, so only particles with high forces move, and the rest are relatively still
		// untill the highest forces get down to their level?

		const float alpha = 8000.f;        // Learning rate (can be tuned)
		const float beta1 = 0.9f;          // Decay rate for first moment
		const float beta2 = 0.999f;        // Decay rate for second moment
		const float epsilon = 1e-8f;

		force *= 1e-8;

		// 2. Update Moment Estimates
		const Float3 firstMoment = adamState->firstMoment * beta1 + force * (1 - beta1);
		const Float3 secondMoment = adamState->secondMoment * beta2 + force * force * (1 - beta2);
		adamState->firstMoment = firstMoment;
		adamState->secondMoment = secondMoment;

		// 4. Compute Bias-Corrected Estimates
		const Float3 firstMomentCorrected = firstMoment / (1 - powf(beta1, step));
		const Float3 secondMomentCorrected = secondMoment / (1 - powf(beta2, step));

		const Float3 deltaPos = (firstMomentCorrected / (secondMomentCorrected.sqrtElementwise() + Float3{ epsilon })) * alpha;
		//if (deltaPos.len() > 5.f) {
		//	force.print('F');
		//	deltaPos.print('D');
		//}
		return pos + Coord{ deltaPos * 1e-8f };
	}

//	__device__ static Coord IntegratePositionEM(const Coord& pos, const Float3& force, const float mass, const float dt, float progress/*step/nSteps*/, const Float3& deltaPosPrev) {
//#ifndef ENABLE_INTEGRATEPOSITION
//		return pos;
//#endif
//		const float alpha = 0.15f;
//		const float stepsize = dt * (progress/alpha * expf(1 - progress/alpha)); // Skewed gaussian
//		
//		const float massPlaceholder = 0.01; // Since we dont use velocities, having a different masses would complicate finding the energy minima. We use this placeholder, so all particles have the same "inertia"
//		const Coord deltaCoord = Coord{ (force * (0.5 / massPlaceholder * stepsize*stepsize)).round() };
//
//		// For the final part of EM we regulate the movement heavily
//		const Float3 deltaPos = deltaCoord.toFloat3();
//		if (progress > 0.95f && deltaPos.len() > deltaPosPrev.len() * 0.9f) {
//			return pos + Coord{ deltaPos * (deltaPosPrev.len() * 0.9f) / (deltaPos.len() + 1e-6) };
//		}
//
//		return pos + deltaCoord;
//	}


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
	__device__ static Float3 ForceActivationFunction(const Float3 force, float scalar=1.f) {

		// Handled inf forces by returning a pseudorandom z force based on global thread index
		if (isinf(force.lenSquared())) {
			return GenerateRandomForce();
		}

		if (isnan(force.lenSquared()))
			force.print('A');

		// 1000 [kJ/mol/nm] is a good target for EM. For EM we will scale the forces below this value * 200
		//const float scaleAbove = 1000.f + 30000.f * (1.f-progress);
		//const float alpha = scaleAbove * LIMA / NANO * KILO; // [1/l N/mol]

		// Apply tanh to the magnitude
		const float alpha = 1000.f * KILO * scalar; // [J/mol/nm]

		// Tanh function, ideal for the 
		const float scaledMagnitude = alpha * tanh(force.len()/alpha);

		//const float scaledMagnitude = force.len() / (1.f + force.len() / (2.f * alpha));
		// Scale the original force vector by the ratio of the new magnitude to the original magnitude
		Float3 scaledForce = force * (scaledMagnitude / (force.len() + 1e-8f)); // Avoid division by zero		

		//printf("Force %f %f %f ScaledF %f %f %f\n", force.x, force.y, force.z, scaledForce.x, scaledForce.y, scaledForce.z);
		return scaledForce;
	}

	// TODO: Clean this up, used args
	__device__ inline void LogCompoundData(const CompoundCompact& compound, int totalParticlesUpperbound, CompoundCoords& compound_coords, 
		const float* potE_sum, const Float3& force, Float3& force_LJ_sol, const SimParams& simparams, SimSignals& simsignals, 
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

	__device__ inline void LogSolventData(const BoxParams& boxparams, const float& potE, const NodeIndex& origo, int id, const Coord& relPos, bool solvent_active, 
		const Float3& force, const Float3& velocity, uint32_t step, float* poteBuffer, Float3* trajBuffer, float* velBuffer, int loggingInterval)
	{
		if (step % loggingInterval != 0) { return; }

		if (solvent_active) {
			const int index = DatabuffersDeviceController::GetLogIndexOfParticle(id, boxparams.n_compounds, step, 
				loggingInterval, boxparams.total_particles_upperbound);

			trajBuffer[index] = LIMAPOSITIONSYSTEM::GetAbsolutePositionNM(origo, relPos);
			poteBuffer[index] = potE;
			velBuffer[index] = velocity.len();
		}
	}

	__device__ constexpr bool isOutsideCutoff(const float dist_sq_reciprocal) {
		if constexpr (HARD_CUTOFF) {
			return dist_sq_reciprocal < cutoffNmSquaredReciprocal_device;	//  1. / (CUTOFF_LM * CUTOFF_LM);
		}
		return false;
	}


	template <typename BoundaryCondition>
	__device__ inline void getCompoundHyperpositionsAsFloat3(const NodeIndex& origo_self, const NodeIndex& queryOrigo, const Float3* const queryRelpositions,
		Float3* const output_buffer, Float3& utility_float3, const int n_particles)
	{
		if (threadIdx.x == 0) {
			const NodeIndex querycompound_hyperorigo = BoundaryCondition::applyHyperpos_Return(origo_self, queryOrigo);
			KernelHelpersWarnings::assertHyperorigoIsValid(querycompound_hyperorigo, origo_self);

			// calc Relative LimaPosition Shift from the origo-shift
			utility_float3 = LIMAPOSITIONSYSTEM_HACK::getRelShiftFromOrigoShift(querycompound_hyperorigo, origo_self).ToRelpos();
		}
		__syncthreads();

		if (threadIdx.x < n_particles) {
			output_buffer[threadIdx.x] = queryRelpositions[threadIdx.x] + utility_float3;
		}
		__syncthreads();
	}

	template <typename BondType, int max_bondtype_in_compound>
	__device__ BondType* LoadBonds(char* utility_buffer, const BondType* const source, int nBondsToLoad) {
		static_assert(cbkernel_utilitybuffer_size >= sizeof(BondType) * max_bondtype_in_compound, "Utilitybuffer not large enough for bondtype");
		BondType* bonds = (BondType*)utility_buffer;

		auto block = cooperative_groups::this_thread_block();
		cooperative_groups::memcpy_async(block, bonds, source, sizeof(BondType) * nBondsToLoad);
		cooperative_groups::wait(block);

		return bonds;
	}

};

