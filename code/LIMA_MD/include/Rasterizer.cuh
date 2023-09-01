#pragma once

#include "cuda_runtime.h"

#include "LimaTypes.cuh"
#include "Simulation.cuh"
//#include "Forcefield.cuh"

const int RAS_THREADS_PER_BLOCK = 64;


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
	
	std::vector<RenderBall> render(Simulation* simulation);

	int solvent_offset = 0;

private:
	/// <summary>	/// Returns a pointer to a list of atoms on the device	/// </summary>
	RenderAtom* getAllAtoms(Simulation* simulation);

	std::vector<RenderBall> processAtoms(RenderAtom* atoms, Simulation* simulation);


	int n_threadblocks = 0;
};
