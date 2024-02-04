#ifndef __linux__


#include "Rasterizer.cuh"
#include "EngineUtils.cuh"
#include "Utilities.h"


#include <algorithm>
#include <cuda/std/cmath>



std::vector<RenderBall> Rasterizer::render(const Float3* positions, const std::vector<Compound>& compounds, const BoxParams& boxparams, int64_t step, Float3 camera_normal) {

    LIMA_UTILS::genericErrorCheck("Error before renderer");
	RenderAtom* atoms_dev = getAllAtoms(positions, compounds, boxparams, step);

    std::vector<RenderBall> balls_host = processAtoms(atoms_dev, boxparams.total_particles_upperbound, boxparams.dims.x);
    cudaFree(atoms_dev);

    std::sort(balls_host.begin(), balls_host.end(), [camera_normal](const RenderBall& a, const RenderBall& b) {
        float dotA = a.pos.x * camera_normal.x + a.pos.y * camera_normal.y + a.pos.z * camera_normal.z;
        float dotB = b.pos.x * camera_normal.x + b.pos.y * camera_normal.y + b.pos.z * camera_normal.z;

        // Compare dot products for sorting
        return !a.disable && dotA > dotB;
        });


    LIMA_UTILS::genericErrorCheck("Error after renderer");

	return balls_host;
}




__global__ void loadCompoundatomsKernel(RenderAtom* atoms, const int step, const Float3* positions, const Compound* compounds);
__global__ void loadSolventatomsKernel(const Float3* positions, int n_compounds, int n_solvents, RenderAtom* atoms);
__global__ void processAtomsKernel(RenderAtom* atoms, RenderBall* balls, float box_len_nm);

RenderAtom* Rasterizer::getAllAtoms(const Float3* positions, const std::vector<Compound>& compounds, const BoxParams& boxparams, int64_t step) {
	RenderAtom* atoms;
	cudaMallocManaged(&atoms, sizeof(RenderAtom) * boxparams.total_particles_upperbound);



  
    Float3* positions_dev;
    cudaMallocManaged(&positions_dev, sizeof(Float3) * boxparams.total_particles_upperbound);
    cudaMemcpy(positions_dev, positions, sizeof(Float3) * boxparams.total_particles_upperbound, cudaMemcpyHostToDevice);

    
    Compound* compounds_dev;
    cudaMallocManaged(&compounds_dev, sizeof(Compound) * boxparams.n_compounds);
    cudaMemcpy(compounds_dev, compounds.data(), sizeof(Compound) * boxparams.n_compounds, cudaMemcpyHostToDevice);


    if (boxparams.n_compounds > 0) {
        loadCompoundatomsKernel << <boxparams.n_compounds, MAX_COMPOUND_PARTICLES >> > (atoms, step, positions_dev, compounds_dev);
    }
	if (boxparams.n_solvents > 0) {
		loadSolventatomsKernel << < boxparams.n_solvents / 128+1, 128 >> > (positions_dev, boxparams.n_compounds, boxparams.n_solvents, atoms);   // TODO: This nsol/128 is wrong, if its not a multiple of 128
	}


    cudaFree(positions_dev);
    cudaFree(compounds_dev);
	return atoms;
}


std::vector<RenderBall> Rasterizer::processAtoms(RenderAtom* atoms, int total_particles_upperbound, float box_len_nm) {
    RenderBall* balls_device;
    cudaMalloc(&balls_device, sizeof(RenderBall) * total_particles_upperbound);
    processAtomsKernel <<< total_particles_upperbound/128+1, 128 >>> (atoms, balls_device, box_len_nm);
    LIMA_UTILS::genericErrorCheck("Error during rendering");

    std::vector<RenderBall> balls_host(total_particles_upperbound);
    cudaMemcpy(balls_host.data(), balls_device, sizeof(RenderBall) * total_particles_upperbound, cudaMemcpyDeviceToHost);

    cudaFree(balls_device);

    return balls_host;
}







__device__ ATOM_TYPE RAS_getTypeFromIndex(int atom_index) {
    switch (atom_index)
    {
    case 0:
        return ATOM_TYPE::SOL;
    case 1:
        return ATOM_TYPE::C;
    case 2:
        return ATOM_TYPE::O;
    case 3:
        return ATOM_TYPE::N;
    case 4:
        return ATOM_TYPE::H;
    case 5: 
        return ATOM_TYPE::P;
    case 6:
        return ATOM_TYPE::SOL;
    default:
        return ATOM_TYPE::NONE;
    }
}

__device__ ATOM_TYPE RAS_getTypeFromMass(double mass) {
    mass *= 1000.f;   //convert to g
    if (mass < 4)
        return ATOM_TYPE::H;
    if (mass < 14)
        return ATOM_TYPE::C;
    if (mass < 15)
        return ATOM_TYPE::N;
    if (mass < 18)
        return ATOM_TYPE::O;
    if (mass < 32)
        return ATOM_TYPE::P;
    return ATOM_TYPE::NONE;
}

__device__ Int3 getColor(ATOM_TYPE atom_type) {
    switch (atom_type)
    {
    case ATOM_TYPE::SOL:
        return Int3(0x03, 0xa9, 0xf4);
    case ATOM_TYPE::H:
        return Int3(0xF1, 0xF1, 0xF1);
    case ATOM_TYPE::O:
        return Int3(0xE0, 0x20, 0x20);
    case ATOM_TYPE::C:
        return Int3(0x30, 0x10, 0x90);
    case ATOM_TYPE::P:
        return Int3(0xFC, 0xF7, 0x5E);
    case ATOM_TYPE::N:
        return Int3(0x2E, 0x8B, 0x57);
    case ATOM_TYPE::NONE:
        return Int3(0xF2, 0xE5, 0xD9);     
    default:
        return Int3(0, 0, 0);
    }
}

//https://en.wikipedia.org/wiki/Van_der_Waals_radius
__device__ float getRadius(ATOM_TYPE atom_type) {
    switch (atom_type)
    {
    case ATOM_TYPE::H:
        return 0.109f * 0.25f;   // Make smaller for visibility
    case ATOM_TYPE::C:
        return 0.17f;
    case ATOM_TYPE::N:
        return 0.155f;
    case ATOM_TYPE::O:
        return 0.152f;
    case ATOM_TYPE::SOL:
        return 0.152f * 0.4f;   // Make smaller for visibility
    case ATOM_TYPE::P:
        return 0.18f;
    case ATOM_TYPE::NONE:
        return 1.f;
    default:
        return 1.f;
    }
}








/// <summary>
/// 
/// </summary>
/// <param name="box"></param>
/// <param name="atoms"></param>
/// <param name="step"></param>
/// <param name="positions">Absolute positions in nm of all particles (compounds first then solvents)</param>
/// <returns></returns>
__global__ void loadCompoundatomsKernel(RenderAtom* atoms, const int step, const Float3* positions, const Compound* compounds) {

    const int local_id = threadIdx.x;
    const int compound_id = blockIdx.x;
    const int global_id = threadIdx.x + blockIdx.x * blockDim.x;

    //Compound* compound = &box->compounds[compound_id];
    const Compound* compound = &compounds[compound_id];
    
    if (local_id < compound->n_particles) {

        RenderAtom atom{};
        atom.pos = positions[global_id];
        atom.mass = SOLVENT_MASS;                                                         // TEMP
        atom.atom_type = RAS_getTypeFromIndex(compound->atom_color_types[local_id]);
        atom.color = getColor(atom.atom_type);
        //ATOM_TYPE typetemp = static_cast<ATOM_TYPE>(compound_id % 7);
        //atom.color = getColor(typetemp);
        atoms[global_id] = atom;
    }
    else {
        atoms[global_id].atom_type = ATOM_TYPE::NONE;
    }
}

__global__ void loadSolventatomsKernel(const Float3* positions, int n_compounds, int n_solvents, RenderAtom* atoms) 
{
    const int solvent_index = blockIdx.x * blockDim.x + threadIdx.x;
    const int particle_index = n_compounds * MAX_COMPOUND_PARTICLES + solvent_index;

    if (solvent_index < n_solvents) {

		RenderAtom atom{};
        atom.pos = positions[particle_index];
		atom.mass = SOLVENT_MASS;
		atom.atom_type = SOL;

		// Debug
		//float velocity = (atom.pos - SolventCoord{ solventblock_prev->origo, solventblock_prev->rel_pos[threadIdx.x] }.getAbsolutePositionLM()).len();
  //      const float velocity = 1.f;
  //      const float point1nm = NANO_TO_LIMA * 0.1f;
		//const float color_scalar = velocity / point1nm * 255.f;
		//const uint8_t color_red = static_cast<uint8_t>(cuda::std::__clamp_to_integral<uint8_t, float>(color_scalar));
		//atom.color = Int3(color_red, 0, 255 - color_red);
        atom.color = Int3(0, 0, 255);

        atoms[particle_index] = atom;
    }
}

__global__ void processAtomsKernel(RenderAtom* atoms, RenderBall* balls, float box_len_nm) { 
    const int index = threadIdx.x + blockIdx.x * blockDim.x;
    
    RenderAtom atom = atoms[index];

    
    //atom.color = getColor(atom.atom_type);

    atom.radius = getRadius(atom.atom_type) * 1.f/box_len_nm;       // [nm]

    // Convert units to normalized units for OpenGL
    atom.pos = atom.pos / box_len_nm - 0.5f;    // normalize from -0.5->0.5
    //atom.pos *= 1.4f;                           // scale a bit up

    //for (int dim = 0; dim < 3; dim++) {
    //    *atom.pos.placeAt(dim) = (atom.pos.at(dim) / BOX_LEN_NM - 0.5f);
    //}
    
    const RenderBall ball(atom.pos, atom.radius, atom.color, atom.atom_type);
    balls[index] = ball;
}
#endif
