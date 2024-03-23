#pragma once

#include<iostream>
//#include <cmath>

#include "LimaTypes.cuh"
#include "Constants.h"
#include "Simulation.cuh"
#include "Forcefield.cuh"
#include "EngineUtilsWarnings.cuh"
#include "BoundaryConditionPublic.h"
#include "LimaPositionSystem.cuh"
#include "LimaTypes.cuh"
#include "Constants.h"
#include "Bodies.cuh"


//#include <cuda.h>
//#include <device_launch_parameters.h>
//#include <cuda_runtime_api.h>

//#include <device_functions.h>

//#include "cuda/std/cmath"
//#include "cuda/std//utility"

namespace CPPD {
	__device__ __host__ constexpr int32_t ceil(float num) {
		return (static_cast<float>(static_cast<int32_t>(num)) == num)
			? static_cast<int32_t>(num)
			: static_cast<int32_t>(num) + ((num > 0) ? 1 : 0);
	}

	template <typename T>
	__device__ __host__ static constexpr T max(const T l, const T r) {
		return r > l ? r : l;
	}

	template <typename T>
	__device__ __host__ static T min(const T l, const T r) {
		return r < l ? r : l;
	}

	__device__ __host__ static int32_t abs(const int32_t val) {
		return val < 0 ? -val : val;
	}
}






namespace LIMALOGSYSTEM {

	// Same as below, but we dont expect to be an even interval
	static constexpr int64_t getMostRecentDataentryIndex(int64_t step, int loggingInterval) {
		return step / loggingInterval;
	}

	static constexpr int64_t getNIndicesBetweenSteps(int64_t from, int64_t to, int loggingInterval) {
		return getMostRecentDataentryIndex(to, loggingInterval) - getMostRecentDataentryIndex(from, loggingInterval);
	}

	__device__ __host__ static constexpr int64_t getDataentryIndex(int64_t step, int loggingInterval) {
		if (step % loggingInterval != 0) {
			//throw std::runtime_error("This step was not expected, there is no equivalent entryindex for it.");	// TODO Maybe then return previous valid entryindex? FIXFIX DANGER
			return 0;
		}
		assert(step % loggingInterval == 0);
		return step / loggingInterval;
	}
};


namespace EngineUtils {



	//// For solvents, compound_id = n_compounds and particle_id = solvent_index
	//__device__ static constexpr int64_t getLoggingIndexOfParticle(uint32_t step, uint32_t total_particles_upperbound, uint32_t compound_id, uint32_t particle_id_local, int loggingInterval) {

	//	const int64_t steps_since_transfer = (step % STEPS_PER_LOGTRANSFER);
	//	//const int64_t step_offset = steps_since_transfer * total_particles_upperbound;
		//const int64_t step_offset = LIMALOGSYSTEM::getDataentryIndex(steps_since_transfer, loggingInterval) * total_particles_upperbound;
	//	const int compound_offset = compound_id * MAX_COMPOUND_PARTICLES;
	//	return step_offset + compound_offset + particle_id_local;
	//}

	template <typename BoundaryCondition>
	__device__ int static getNewBlockId(const NodeIndex& transfer_direction, const NodeIndex& origo) {
		NodeIndex new_nodeindex = transfer_direction + origo;
		BoundaryCondition::applyBC(new_nodeindex);
		return SolventBlocksCircularQueue::get1dIndex(new_nodeindex);
	}

	//__device__ static SolventBlockTransfermodule* getTransfermoduleTargetPtr(SolventBlockTransfermodule* transfermodule_array, int blockId, const NodeIndex& transfer_direction) {
	//	NodeIndex new_nodeindex = SolventBlocksCircularQueue::get3dIndex(blockId) + transfer_direction;
	//	LIMAPOSITIONSYSTEM::applyBC(new_nodeindex);
	//	const int index = SolventBlocksCircularQueue::get1dIndex(new_nodeindex);
	//	return &transfermodule_array[index];
	//}


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

	__device__ inline void LogCompoundData(const CompoundCompact& compound, Box* box, CompoundCoords& compound_coords, 
		float* potE_sum, const Float3& force, Float3& force_LJ_sol, const SimParams& simparams, SimSignals& simsignals, DatabuffersDevice* databuffers, const float speed)
	{
		if (threadIdx.x >= compound.n_particles) { return; }

		if (simsignals.step % simparams.data_logging_interval != 0) { return; }

		//const int64_t index = EngineUtils::getLoggingIndexOfParticle(simsignals.step, box->boxparams.total_particles_upperbound, blockIdx.x, threadIdx.x, simparams.data_logging_interval);
		const int index = databuffers->GetLogIndexOfParticle(threadIdx.x, blockIdx.x, simsignals.step, simparams.data_logging_interval);
		databuffers->traj_buffer[index] = LIMAPOSITIONSYSTEM::getAbsolutePositionNM(compound_coords.origo, compound_coords.rel_positions[threadIdx.x]); //LIMAPOSITIONSYSTEM::getGlobalPosition(compound_coords);
		databuffers->potE_buffer[index] = *potE_sum;
		databuffers->vel_buffer[index] = speed;

		EngineUtilsWarnings::logcompoundVerifyVelocity(compound, simparams, simsignals, compound_coords, force, speed);
	}

	__device__ inline void LogSolventData(Box* box, const float& potE, const SolventBlock& solventblock, bool solvent_active, 
		const Float3& force, const Float3& velocity, uint32_t step, DatabuffersDevice* databuffers, int loggingInterval)
	{
		if (step % loggingInterval != 0) { return; }


		if (solvent_active) {
			//const int64_t index = EngineUtils::getLoggingIndexOfParticle(step, box->boxparams.total_particles_upperbound, box->boxparams.n_compounds, solventblock.ids[threadIdx.x], loggingInterval);
			const int index = databuffers->GetLogIndexOfParticle(solventblock.ids[threadIdx.x], box->boxparams.n_compounds, step, loggingInterval);


			databuffers->traj_buffer[index] = LIMAPOSITIONSYSTEM::getAbsolutePositionNM(solventblock.origo, solventblock.rel_pos[threadIdx.x]);
			databuffers->potE_buffer[index] = potE;
			databuffers->vel_buffer[index] = velocity.len();
			

#ifdef USEDEBUGF3
			const auto debug_index = (box->step * box->total_particles_upperbound + box->n_compounds * MAX_COMPOUND_PARTICLES + solventblock.ids[threadIdx.x]) * DEBUGDATAF3_NVARS;
			//box->debugdataf3[debug_index] = Float3(solventblock.ids[threadIdx.x] * 10 + 1.f, solventblock.ids[threadIdx.x] * 10 + 2.f, solventblock.ids[threadIdx.x] * 10 + 3.f);
			box->debugdataf3[debug_index] = force;
			box->debugdataf3[debug_index + 1] = velocity;
			box->debugdataf3[debug_index + 2] = SolventBlockHelpers::extractAbsolutePositionLM(solventblock) / NANO_TO_LIMA;
#endif
		}
	}



	//// Slow function, never for device
	//__host__ static float calcDistance(const NodeIndex& o1, const Coord& relpos1, const NodeIndex& o2, const Coord& relpos2) {
	//	auto pos1 = LIMAPOSITIONSYSTEM::getAbsolutePositionNM(o1, relpos1);
	//	auto pos2 = LIMAPOSITIONSYSTEM::getAbsolutePositionNM(o2, relpos2);
	//	return (pos1 - pos2).len();
	//}

	__device__ static bool constexpr isOutsideCutoff(const float dist_sq_reciprocal) {
		if constexpr (HARD_CUTOFF) {
			constexpr float threshold = 1. / (CUTOFF_LM * CUTOFF_LM);
			if (dist_sq_reciprocal < threshold) {
				return true;
			}
		}
		return false;
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
	void findAllNearestSolventSolvent(SolventBlocksCircularQueue* queue, size_t n_solvents, std::vector<float>& out);
}

// LIMA algorithm Library
namespace LAL {
	//__device__ constexpr int getBlellochTablesize(int n) {
	//	const float nf = static_cast<float>(n);
	//	return CPPD::ceil(nf * std::log2f(nf) * 2.f);
	//}
	// 



}
