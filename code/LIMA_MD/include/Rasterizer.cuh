#pragma once

#include "cuda_runtime.h"

#include "LimaTypes.cuh"
#include "Simulation.cuh"
#include <limits>


static_assert(sizeof(RenderAtom) == 32, "RenderAtom size is not 32 bytes, and thus risc being packed wrongly for GLSL");


class Rasterizer {
public:
	Rasterizer() {};
	~Rasterizer();

	void render(const Float3* positions,
		const std::vector<Compound>& compounds, const BoxParams&, int64_t step, Float3 camera_normal, ColoringMethod coloringmethod, RenderAtom* renderAtoms);


private:
	
	void getAllAtoms(const Float3* positions, const std::vector<Compound>& compounds, 
		const BoxParams& boxparams, int64_t step, ColoringMethod coloringMethod, RenderAtom* renderAtoms);


	//RenderAtom* atoms_dev = nullptr;
	Float3* positions_dev = nullptr;
	Compound* compounds_dev = nullptr;

	bool isInitialized = false;
	void initialize(const BoxParams& boxparams, const std::vector<Compound>& compounds);

};
