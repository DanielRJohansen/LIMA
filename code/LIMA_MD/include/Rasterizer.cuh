#pragma once

#include "LimaTypes.cuh"
#include "cuda_runtime.h"
#include "Simulation.cuh"
//#include "Forcefield.cuh"

const int RAS_THREADS_PER_BLOCK = 64;


struct RenderAtom {
	Float3 pos;
	float mass = 0;
	float radius = 0;
	Int3 color{};
	ATOM_TYPE atom_type = ATOM_TYPE::NONE;
};



class Rasterizer {
public:
	Rasterizer() {};
	
	RenderBall* render(Simulation* simulation);

	//int actual_n_particles;
	//int n_particles_upperbound;
	int solvent_offset = 0;

private:
	RenderAtom* getAllAtoms(Simulation* simulation);
	void sortAtoms(RenderAtom* atoms, int dim);
	RenderBall* processAtoms(RenderAtom* atoms, Simulation* simulation);


	int n_threadblocks = 0;
	//int actual_n_particles;
};
