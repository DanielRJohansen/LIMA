#pragma once

#include "LimaTypes.cuh"
#include "EngineBodies.cuh"
#include "DeviceAlgorithms.cuh"

#include "KernelConstants.cuh"

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
			const float box_padding = 0.5f;	// The dist to the box edges (from compound center) we want to enforce, so switching to PBC wont cause immediate collisions

			const float boxlenHalfNM = boxSize_device.boxSizeNM_f / 2.f;
			const float dist_x = LAL::max(std::abs(boxlenHalfNM - avg_compound_position_nm.x) - boxlenHalfNM + box_padding, 0.f);
			const float dist_y = LAL::max(std::abs(boxlenHalfNM - avg_compound_position_nm.y) - boxlenHalfNM + box_padding, 0.f);


			// Constant force
			float force_x{}, force_y{};

			const float const_factor = 0.00001f;
			force_x += const_factor;
			force_y += const_factor;

			// Linear force
			const float lin_factor = 0.00002f;
			force_x += dist_x * lin_factor;
			force_y += dist_y * lin_factor;

			// Exp force
			/*const float exp_factor = 0.000000;
			force_x += dist_x * dist_x * exp_factor;
			force_y += dist_y * dist_y * exp_factor;*/

			// Apply signed force
			if (avg_compound_position_nm.x > boxlenHalfNM)
				force_x *= -1.f;
			if (avg_compound_position_nm.y > boxlenHalfNM)
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

			//if (std::abs(avg_compound_force_z) > 10.f)
			//	printf("force %f \n", avg_compound_force_z);

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

		utilitybuffer_f[threadIdx.x] = 0;
		utilitybuffer_f3[threadIdx.x] = Float3{ 0 };
		__syncthreads();

		// Compute the avg position
		{
			utilitybuffer_f3[threadIdx.x] = compound_positions[threadIdx.x] * LIMA_TO_NANO;
			__syncthreads();
			distributedSummation(utilitybuffer_f3, n_particles);	// Get summed relpos in nm at pos[0]
			if (threadIdx.x == 0) {
				utilitybuffer_f3[0] *= 1.f / static_cast<float>(n_particles);	// Get average relpos in nm
				*avg_abspos_nm = LIMAPOSITIONSYSTEM::GetAbsolutePositionNM(compound_origo, Coord{ utilitybuffer_f3[0] });
			}
		}
		__syncthreads();
		// Compute the avg force
		{
			utilitybuffer_f[threadIdx.x] = force.z;
			__syncthreads();
			distributedSummation(utilitybuffer_f, n_particles);
			if (threadIdx.x == 0) {
				*avg_force_z = utilitybuffer_f[0] / static_cast<float>(n_particles);
			}
		}
		__syncthreads();

		if (threadIdx.x < n_particles)
			SupernaturalForces::_applyHorizontalSqueeze(*avg_abspos_nm, *avg_force_z, force, mass);
		__syncthreads();
		*avg_abspos_nm = Float3{ 0 };
		*avg_force_z = 0;
		utilitybuffer_f[threadIdx.x] = 0;
		utilitybuffer_f3[threadIdx.x] = Float3{ 0 };
		__syncthreads();
	}


	// Overwrites the force that is given as an argument to the function
	__global__ void ApplyHorizontalSqueeze(SimulationDevice* simDev) {
		__shared__ Float3 avg_abspos_nm;
		__shared__ float avg_force_z;
		__shared__ Float3 relPosNm[MAX_COMPOUND_PARTICLES];
		__shared__ float forcesZ[MAX_COMPOUND_PARTICLES];
		__shared__ Float3 origo;

		relPosNm[threadIdx.x] = Float3{ 0 };
		forcesZ[threadIdx.x] = 0;

		int nParticles = simDev->box->compounds[blockIdx.x].n_particles;

		__syncthreads();


		// Compute the avg position
		{
			const auto const coords = simDev->box->compoundcoordsCircularQueue->getCoordarrayRef(simDev->signals->step, blockIdx.x);
			LIMAPOSITIONSYSTEM::LoadCompoundPositionAsLm(coords, origo, relPosNm, nParticles);
			__syncthreads();
			relPosNm[threadIdx.x] = relPosNm[threadIdx.x] * LIMA_TO_NANO + origo;
			__syncthreads();

			distributedSummation(relPosNm, MAX_COMPOUND_PARTICLES);	// Get summed relpos in nm at pos[0]
			if (threadIdx.x == 0) {
				avg_abspos_nm = relPosNm[0] / static_cast<float>(nParticles);	// Get average relpos in nm
			}
		}

		// Compute the avg force
		{
			if (threadIdx.x < nParticles) {
				forcesZ[threadIdx.x] = simDev->box->compounds[blockIdx.x].forces_interim[threadIdx.x].z;
			}
			__syncthreads();				

			distributedSummation(forcesZ, nParticles);

			if (threadIdx.x == 0) {
				avg_force_z = 0.f;// forcesZ[0] / static_cast<float>(nParticles);
			}
		}
		__syncthreads();

		
		if (threadIdx.x < nParticles) {
			Float3 force{};
			const float mass = forcefield_device.particle_parameters[simDev->box->compounds[blockIdx.x].atom_types[threadIdx.x]].mass;
			SupernaturalForces::_applyHorizontalSqueeze(avg_abspos_nm, avg_force_z, force, mass);
			if (force.len() > 10.f ) {
				force.print('f');
				avg_abspos_nm.print('p');
			}

			simDev->box->compounds[blockIdx.x].forces_interim[threadIdx.x] += force;
			//if (simDev->box->compounds[blockIdx.x].forces_interim[threadIdx.x].len() > 10.f) {
			//	simDev->box->compounds[blockIdx.x].forces_interim[threadIdx.x].print('T');
			//}	
		}		
	}
	//__device__ void applyHorizontalChargefield(Float3 posNM, Float3& force, float particleCharge) {
	//	PeriodicBoundaryCondition::applyBCNM(posNM);	// TODO: Use generic BC
	//	const float distFromMidPlane = posNM.x - (boxSize_device.boxSizeNM_f / 2.f);

	//	const float dir = distFromMidPlane / std::abs(distFromMidPlane);
	//	const float forceApplied = particleCharge * dir * KILO * KILO * 1000.f;
	//	printf("Force %f\n", forceApplied);
	//	force.x += forceApplied;

	//	if (distFromMidPlane > 5.f)
	//		printf("dist %f charge %f dir %f force %f\n", distFromMidPlane, particleCharge, dir, forceApplied);

	//}
}