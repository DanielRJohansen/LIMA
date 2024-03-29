//#ifndef __linux__


#include "Rasterizer.cuh"
#include "Utilities.h"


#include <algorithm>
#include <cuda/std/cmath>

__global__ void loadCompoundatomsKernel(RenderAtom* atoms, const int step, const Float3* positions, const Compound* compounds, ColoringMethod coloringMethod, float boxLenNM);
__global__ void loadSolventatomsKernel(const Float3* positions, int n_compounds, int n_solvents, RenderAtom* atoms, float boxLenNM);


void Rasterizer::initialize(const BoxParams& boxparams) {
    cudaMallocManaged(&atoms_dev, sizeof(RenderAtom) * boxparams.total_particles_upperbound);
    cudaMallocManaged(&positions_dev, sizeof(Float3) * boxparams.total_particles_upperbound);
    cudaMallocManaged(&compounds_dev, sizeof(Compound) * boxparams.n_compounds);
    renderAtomsHost.resize(boxparams.total_particles_upperbound);

	isInitialized = true;
}

Rasterizer::~Rasterizer() {
    if (isInitialized) {
        cudaFree(positions_dev);
        cudaFree(compounds_dev);
        cudaFree(atoms_dev);
    }
}

/// <summary>	/// Returns a pointer to a list of atoms on the device	/// </summary>
void Rasterizer::getAllAtoms(const Float3* positions, const std::vector<Compound>& compounds, const BoxParams& boxparams, int64_t step, ColoringMethod coloringMethod) {
    cudaMemcpy(positions_dev, positions, sizeof(Float3) * boxparams.total_particles_upperbound, cudaMemcpyHostToDevice);
    cudaMemcpy(compounds_dev, compounds.data(), sizeof(Compound) * boxparams.n_compounds, cudaMemcpyHostToDevice);

    if (boxparams.n_compounds > 0) {
        loadCompoundatomsKernel << <boxparams.n_compounds, MAX_COMPOUND_PARTICLES >> > (atoms_dev, step, positions_dev, compounds_dev, coloringMethod, boxparams.dims.x);
    }
    if (boxparams.n_solvents > 0) {
        loadSolventatomsKernel << < boxparams.n_solvents / 128 + 1, 128 >> > (positions_dev, boxparams.n_compounds, boxparams.n_solvents, atoms_dev, boxparams.dims.x);   // TODO: This nsol/128 is wrong, if its not a multiple of 128
    }
}

const std::vector<RenderAtom>& Rasterizer::render(const Float3* positions, const std::vector<Compound>& compounds, const BoxParams& boxparams, int64_t step, Float3 camera_normal, ColoringMethod coloringMethod) {

    if (!isInitialized) { initialize(boxparams); }

    LIMA_UTILS::genericErrorCheck("Error before renderer");
	getAllAtoms(positions, compounds, boxparams, step, coloringMethod);

    cudaMemcpy(renderAtomsHost.data(), atoms_dev, sizeof(RenderAtom) * boxparams.total_particles_upperbound, cudaMemcpyDeviceToHost);

    std::sort(renderAtomsHost.begin(), renderAtomsHost.end(), [camera_normal](const RenderAtom& a, const RenderAtom& b) {
        float dotA = a.pos.x * camera_normal.x + a.pos.y * camera_normal.y + a.pos.z * camera_normal.z;
        float dotB = b.pos.x * camera_normal.x + b.pos.y * camera_normal.y + b.pos.z * camera_normal.z;

        // Compare dot products for sorting
        return !a.Disabled() && dotA > dotB;
        });

    LIMA_UTILS::genericErrorCheck("Error after renderer\n");

	return renderAtomsHost;
}


__device__ ATOM_TYPE RAS_getTypeFromAtomletter(char atom) {
    switch (atom)
    {
    case 'C':
        return ATOM_TYPE::C;
    case 'O':
        return ATOM_TYPE::O;
    case 'N':
        return ATOM_TYPE::N;
    case 'H':
        return ATOM_TYPE::H;
    case 'P':
        return ATOM_TYPE::P;
    case 'l':
		return ATOM_TYPE::LIMA_CUSTOM;
    default:
        return ATOM_TYPE::NONE;
    }
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
    case ATOM_TYPE::LIMA_CUSTOM:
	    return 0.1f;
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
__global__ void loadCompoundatomsKernel(RenderAtom* atoms, const int step, const Float3* positions, const Compound* compounds, ColoringMethod coloringMethod, float boxLenNM) {

    const int local_id = threadIdx.x;
    const int compound_id = blockIdx.x;
    const int global_id = threadIdx.x + blockIdx.x * blockDim.x;

    //Compound* compound = &box->compounds[compound_id];
    const Compound* compound = &compounds[compound_id];
    
    if (local_id < compound->n_particles) {

        RenderAtom atom{};
        atom.pos = positions[global_id];
        atom.pos = atom.pos / boxLenNM - 0.5f;    // normalize from -0.5->0.5

        //atom.mass = SOLVENT_MASS;                                                         // TEMP
        atom.atom_type = RAS_getTypeFromAtomletter(compound->atomLetters[local_id]);

        atom.radius = getRadius(atom.atom_type) / boxLenNM;

        if (coloringMethod == ColoringMethod::Atomname)
            atom.color = getColor(atom.atom_type);
        else if (coloringMethod == ColoringMethod::Charge) {
            const float chargeNormalized = (compound->atom_charges[local_id] + elementaryChargeToCoulombPerMole)  / (elementaryChargeToCoulombPerMole*2.f);
            atom.color = Int3(255.f * chargeNormalized, 0, 255.f * (1.f - chargeNormalized));
        }

        atoms[global_id] = atom;
    }
    else {
        atoms[global_id].atom_type = ATOM_TYPE::NONE;
    }
}

__global__ void loadSolventatomsKernel(const Float3* positions, int n_compounds, int n_solvents, RenderAtom* atoms, float boxLenNM)
{
    const int solvent_index = blockIdx.x * blockDim.x + threadIdx.x;
    const int particle_index = n_compounds * MAX_COMPOUND_PARTICLES + solvent_index;

    if (solvent_index < n_solvents) {

		RenderAtom atom{};
        atom.pos = positions[particle_index];
        atom.pos = atom.pos / boxLenNM - 0.5f;    // normalize from -0.5->0.5

		atom.mass = SOLVENT_MASS;
		atom.atom_type = SOL;
        atom.radius = getRadius(atom.atom_type) / boxLenNM;
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
