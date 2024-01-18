#pragma once

#include "cuda_runtime.h"

#include "LimaTypes.cuh"
#include "Simulation.cuh"



struct RenderAtom {
	Float3 pos{};		// [nm]
	float mass = 0;		// [kg/mol]
	float radius = 0;	// [??]
	Int3 color{};
	ATOM_TYPE atom_type = ATOM_TYPE::NONE;
};



class Rasterizer {
public:
	Rasterizer() {};
	
	std::vector<RenderBall> render(const Float3* positions,
		const std::vector<Compound>& compounds, const BoxParams&, int64_t step, Float3 camera_normal);


private:
	/// <summary>	/// Returns a pointer to a list of atoms on the device	/// </summary>
	RenderAtom* getAllAtoms(const Float3* positions, 
		const std::vector<Compound>& compounds, const BoxParams& boxparams, int64_t step);

	std::vector<RenderBall> processAtoms(RenderAtom* atoms, int total_particles_upperbound, float boxlen_nm);


};
