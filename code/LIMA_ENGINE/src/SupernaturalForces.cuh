#pragma once

#include "LimaTypes.cuh"
#include "EngineBodies.cuh"

namespace SupernaturalForces {

	// TODO: Move this to LIMA_BASE
	template<typename T>
	void __device__ distributedSummation(T* arrayptr, int array_len) {				// Places the result at pos 0 of input_array
		T temp;			// This is a lazy soluation, but maybe it is also fast? Definitely simple..
		for (int i = 1; i < array_len; i *= 2) {	// Distributed averaging							// Make a generic and SAFER function for this, PLEASE OK??
			if ((threadIdx.x + i) < array_len) {
				temp = arrayptr[threadIdx.x] + arrayptr[threadIdx.x + i];
			}
			__syncthreads();
			arrayptr[threadIdx.x] = temp;
			__syncthreads();
		}
	}



	namespace {// anon namespace
		__device__ void _applyHorizontalSqueeze(const Float3& avg_compound_position_nm, const float& avg_compound_force_z, Float3& particle_force) {
			const float box_padding = 0.3;	// The dist to the box edges (from compound center) we want to enforce, so switching to PBC wont cause immediate collisions

			const float dist_x = BOX_LEN_HALF_NM - avg_compound_position_nm.x;
			const float dist_y = BOX_LEN_HALF_NM - avg_compound_position_nm.y;

			const float factor = 0.00001f;

			const float force_x = std::abs(dist_x) > BOX_LEN_HALF_NM - box_padding
				? dist_x * dist_x * factor * (dist_x < 0 * -1.f) 
				: 0.f;
			const float force_y = std::abs(dist_y) > BOX_LEN_HALF_NM - box_padding
				? dist_y * dist_y * factor * (dist_y < 0 * -1.f)
				: 0.f;

			// Apply squeeze force
			//particle_force.x += force_x;
			//particle_force.y += force_y;
			//if (avg_compound_position_nm.y > BOX_LEN_NM)
			//	particle_force.y -= 0.000001f;


			// Elimitenate z translation
			particle_force.z -= avg_compound_force_z;
		};
	}

	// Overwrites the force that is given as an argument to the function
	__device__ void applyHorizontalSqueeze(Float3* utilitybuffer_f3, float* utilitybuffer_f, char* utilitybuffer, const Float3* const compound_positions, int n_particles, const NodeIndex& compound_origo, Float3& force) {
		static_assert(cbkernel_utilitybuffer_size >= sizeof(Float3) * 2);
		Float3* avg_abspos_nm = reinterpret_cast<Float3*>(utilitybuffer);
		float* avg_force_z = reinterpret_cast<float*>(&utilitybuffer[sizeof(Float3)]);

		// Compute the avg position
		{
			utilitybuffer_f3[threadIdx.x] = compound_positions[threadIdx.x] * LIMA_TO_NANO;
			distributedSummation(utilitybuffer_f3, n_particles);	// Get summed relpos in nm at pos[0]
			if (threadIdx.x == 0) {
				utilitybuffer_f3[0] *= 1.f / static_cast<float>(n_particles);	// Get average relpos in nm
				*avg_abspos_nm = LIMAPOSITIONSYSTEM::getAbsolutePositionNM(compound_origo, Coord{ utilitybuffer_f3[0] });
			}
		}
		__syncthreads();
		// Compute the avg force
		{
			utilitybuffer_f[threadIdx.x] = force.z;
			distributedSummation(utilitybuffer_f, n_particles);
			if (threadIdx.x == 0) {
				*avg_force_z = utilitybuffer_f[0] / static_cast<float>(n_particles);
			}
		}
		__syncthreads();

		SupernaturalForces::_applyHorizontalSqueeze(*avg_abspos_nm, *avg_force_z, force);
	}
}