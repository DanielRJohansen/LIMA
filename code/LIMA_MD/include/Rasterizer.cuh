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

	bool Disabled() const { return atom_type == ATOM_TYPE::NONE; }
};



class Rasterizer {
public:
	Rasterizer() {};
	~Rasterizer();

	const std::vector<RenderAtom>& render(const Float3* positions,
		const std::vector<Compound>& compounds, const BoxParams&, int64_t step, Float3 camera_normal, ColoringMethod coloringmethod);


private:
	
	void getAllAtoms(const Float3* positions, const std::vector<Compound>& compounds, 
		const BoxParams& boxparams, int64_t step, ColoringMethod coloringMethod);


	RenderAtom* atoms_dev = nullptr;
	Float3* positions_dev = nullptr;
	Compound* compounds_dev = nullptr;
	std::vector<RenderAtom> renderAtomsHost;

	bool isInitialized = false;
	void initialize(const BoxParams& boxparams);

};
