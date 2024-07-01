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

#include "KernelWarnings.cuh"
#include "KernelConstants.cuh"

#include <cooperative_groups.h>
#include <cooperative_groups/memcpy_async.h>


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

	__device__ static Float3 SlowHighEnergyParticle(const Float3& velocityNow, const float dt, const float mass, float thermostatScalar) {
		const float kineticEnergy = PhysicsUtils::calcKineticEnergy(velocityNow.len(), mass);
		const float temperature = PhysicsUtils::kineticEnergyToTemperature(kineticEnergy);

		const float temperatureThreshold = 400.f; // [K]

		if (temperature > temperatureThreshold) {
			return velocityNow * 0.9995f;
		}
		else {
			return velocityNow * thermostatScalar;
		}

	}

	template <bool energyMinimize>
	__device__ Float3 ScaleVelocity(const Float3& velocity, const float thermostatScalar, float dt, float mass) {
		if constexpr (energyMinimize)
			return EngineUtils::SlowHighEnergyParticle(velocity, dt, mass, thermostatScalar);
		else
			return velocity * thermostatScalar;
	}




	__device__ inline void LogCompoundData(const CompoundCompact& compound, Box* box, CompoundCoords& compound_coords, 
		float* potE_sum, const Float3& force, Float3& force_LJ_sol, const SimParams& simparams, SimSignals& simsignals, DatabuffersDevice* databuffers, const float speed)
	{
		if (threadIdx.x >= compound.n_particles) { return; }

		if (simsignals.step % simparams.data_logging_interval != 0) { return; }

		const int index = databuffers->GetLogIndexOfParticle(threadIdx.x, blockIdx.x, simsignals.step, simparams.data_logging_interval);
		databuffers->traj_buffer[index] = LIMAPOSITIONSYSTEM::getAbsolutePositionNM(compound_coords.origo, compound_coords.rel_positions[threadIdx.x]); 
		databuffers->potE_buffer[index] = *potE_sum;
		databuffers->vel_buffer[index] = speed;

		EngineUtilsWarnings::logcompoundVerifyVelocity(compound, simparams, simsignals, compound_coords, force, speed);
	}

	__device__ inline void LogSolventData(Box* box, const float& potE, const SolventBlock& solventblock, bool solvent_active, 
		const Float3& force, const Float3& velocity, uint32_t step, DatabuffersDevice* databuffers, int loggingInterval)
	{
		if (step % loggingInterval != 0) { return; }


		if (solvent_active) {
			const int index = databuffers->GetLogIndexOfParticle(solventblock.ids[threadIdx.x], box->boxparams.n_compounds, step, loggingInterval);

			databuffers->traj_buffer[index] = LIMAPOSITIONSYSTEM::getAbsolutePositionNM(solventblock.origo, solventblock.rel_pos[threadIdx.x]);
			databuffers->potE_buffer[index] = potE;
			databuffers->vel_buffer[index] = velocity.len();
		}
	}

	__device__ static bool constexpr isOutsideCutoff(const float dist_sq_reciprocal) {
		if constexpr (HARD_CUTOFF) {
			constexpr float threshold = 1. / (CUTOFF_LM * CUTOFF_LM);
			if (dist_sq_reciprocal < threshold) {
				return true;
			}
		}
		return false;
	}










	template <typename BoundaryCondition>
	__device__ inline void getCompoundHyperpositionsAsFloat3(const NodeIndex& origo_self, const CompoundCoords* const querycompound,
		void* output_buffer, Float3& utility_float3, const int n_particles)
	{
		if (threadIdx.x == 0) {
			const NodeIndex querycompound_hyperorigo = LIMAPOSITIONSYSTEM::getHyperNodeIndex<BoundaryCondition>(origo_self, querycompound->origo);
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
			const Coord queryparticle_coord = ((Coord*)output_buffer)[threadIdx.x];
			((Float3*)output_buffer)[threadIdx.x] = queryparticle_coord.toFloat3() + relshift;
		}
	}

	

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

