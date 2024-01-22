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
		__device__ void _applyHorizontalSqueeze(const Float3& avg_compound_position_nm, const float& avg_compound_force_z, Float3& particle_force, float particle_mass) {
			const float box_padding = 1.4f;	// The dist to the box edges (from compound center) we want to enforce, so switching to PBC wont cause immediate collisions

			const float dist_x = CPPD::max(std::abs(BOX_LEN_HALF_NM - avg_compound_position_nm.x) - BOX_LEN_HALF_NM + box_padding, 0.f);
			const float dist_y = CPPD::max(std::abs(BOX_LEN_HALF_NM - avg_compound_position_nm.y) - BOX_LEN_HALF_NM + box_padding, 0.f);


			// Constant force
			float force_x{}, force_y{};

			const float const_factor = 0.00001f;
			force_x += const_factor;
			force_y += const_factor;

			// Linear force
			const float lin_factor = 0.00003f;
			force_x += dist_x * lin_factor;
			force_y += dist_y * lin_factor;

			// Exp force
			const float exp_factor = 0.000001;
			force_x += dist_x * dist_x * exp_factor;
			force_y += dist_y * dist_y * exp_factor;

			// Apply signed force
			if (avg_compound_position_nm.x > BOX_LEN_HALF_NM)
				force_x *= -1.f;
			if (avg_compound_position_nm.y > BOX_LEN_HALF_NM)
				force_y *= -1.f;

			float mass_factor = particle_mass * 100.f;	// If we scale the forces by their inverted mass, we avoid the problem of lighter molecules being pushed faster than the heavier
			if (mass_factor < 0.1f || mass_factor > 4.f) {
				printf("mass factor %f mass %f\n", mass_factor, particle_mass);
			}	
			//mass_factor = 1.f;

			// Apply squeeze force
			particle_force.x += force_x * mass_factor;
			particle_force.y += force_y * mass_factor;
			//if (avg_compound_position_nm.y > BOX_LEN_NM)
			//	particle_force.y -= 0.000001f;


			// Elimitenate z translation
			particle_force.z -= avg_compound_force_z;// *mass_factor;
		};
	}

	// Overwrites the force that is given as an argument to the function
	__device__ void applyHorizontalSqueeze(Float3* utilitybuffer_f3, float* utilitybuffer_f, char* utilitybuffer, const Float3* const compound_positions, 
		int n_particles, const NodeIndex& compound_origo, Float3& force, float mass) {
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

		if (threadIdx.x < n_particles)
			SupernaturalForces::_applyHorizontalSqueeze(*avg_abspos_nm, *avg_force_z, force, mass);
		__syncthreads();
	}
}