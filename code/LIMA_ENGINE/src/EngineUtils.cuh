#pragma once


#include "LimaTypes.cuh"
#include "Constants.h"
#include "Simulation.cuh"
#include "SimulationData.h"
#include "EngineUtilsWarnings.cuh"
#include "LimaPositionSystem.cuh"
#include "LimaTypes.cuh"
#include "Constants.h"
#include "Bodies.cuh"
#include "BoxGrid.cuh"
#include "KernelWarnings.cuh"

#include <cooperative_groups.h>
#include <cooperative_groups/memcpy_async.h>

namespace EngineUtils {

	template <typename BoundaryCondition>
	__device__ int static getNewBlockId(const NodeIndex& transferDirection, const NodeIndex& origo, const Int3& boxSize) {
		NodeIndex newNodeIndex = transferDirection + origo;
		BoundaryCondition::applyBC(newNodeIndex, boxSize);
		return BoxGrid::Get1dIndex(newNodeIndex, boxSize);
	}

	// returns pos_tadd1
	__device__ static Coord integratePositionVVS(const Coord& pos, const Float3& vel, const Float3& force, const float mass, const float dt) {
		if constexpr (!ENABLE_INTEGRATEPOSITION) {
			return pos;
		}

		const Coord pos_tadd1 = pos + Coord{ (vel * dt + force * (0.5f / mass * dt * dt)) };				// precise version
		return pos_tadd1;
	}
	constexpr static Float3 IntegratePositionVVS(const Float3& pos, const Float3& vel, const Float3& force, const float mass, const float dt) {
		if constexpr (!ENABLE_INTEGRATEPOSITION) {
			return pos;
		}

		const Float3 pos_tadd1 = pos + (vel * dt + force * (0.5f / mass * dt * dt));				// precise version
		return pos_tadd1;
	}
	__device__ static Float3 integrateVelocityVVS(const Float3& vel_tsub1, const Float3& force_tsub1, const Float3& force, const float dt, const float mass) {
		const Float3 vel = vel_tsub1 + (force + force_tsub1) * (dt * 0.5f / mass);
		return vel;
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
	__device__ inline Float3 GenerateRandomForce(int pidGlobal) {
		unsigned int seed = pidGlobal;

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
	__device__ static Float3 ForceActivationFunction(int pidGlobal /*Used as a random-seed*/, const Float3 force, float scalar = 1.f) {

		// Handled inf forces by returning a pseudorandom z force based on global thread index
		if (isinf(force.lenSquared())) {
			return GenerateRandomForce(pidGlobal);
		}

		if (isnan(force.lenSquared())) {
			force.print('A');
			asm("trap;"); // Force the kernel to crash
		}

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

	__device__ inline void LogPclusterData(int pcId, int pidInPclusters, int64_t step, int data_logging_interval, Float3 position, float potential, Float3 force, float speed, int totalParticlesUpperbound,
		Float3* trajBuffer, float* potEBuffer, float* velocityBuffer, Float3* forceBuffer) {
		//if (threadIdx.x >= compound.n_particles) { return; }

		if (data_logging_interval == 0 || step % data_logging_interval != 0) { return; }

		const size_t index = DatabuffersDeviceController::GetLogIndexOfParticle(pidInPclusters, pcId, step, data_logging_interval, totalParticlesUpperbound);
		trajBuffer[index] = position;
		potEBuffer[index] = potential;
		velocityBuffer[index] = speed;
		forceBuffer[index] = force;
	}


	__device__ constexpr bool isOutsideCutoff(const float distSq, const float cutoffNmSquared) {
		if constexpr (HARD_CUTOFF) {
			return distSq > cutoffNmSquared;	// (CUTOFF_LM * CUTOFF_LM);
		}
		return false;
	}

    __device__ constexpr bool isOutsideCutoff_recip(const float distSqReciprocal, const float cutoffNmSquaredReciprocal) {
        if constexpr (HARD_CUTOFF) {
            return distSqReciprocal < cutoffNmSquaredReciprocal;
        }
        return false;
    }

};

