//#ifndef __linux__


#include "Rasterizer.cuh"
#include "Utilities.h"


#include <algorithm>
#include <cuda/std/cmath>

__global__ void loadCompoundatomsKernel(RenderAtom* atoms, const int step, const Float3* positions, const Compound* compounds, ColoringMethod coloringMethod, float boxLenNM);
__global__ void loadSolventatomsKernel(const Float3* positions, int n_compounds, int n_solvents, RenderAtom* atoms, float boxLenNM);


void Rasterizer::initialize(const BoxParams& boxparams) {
    cudaMallocManaged(&atoms_dev, sizeof(RenderAtom) * boxparams.total_particles);
    cudaMallocManaged(&positions_dev, sizeof(Float3) * boxparams.total_particles_upperbound);
    cudaMallocManaged(&compounds_dev, sizeof(Compound) * boxparams.n_compounds);
    renderAtomsHost.resize(boxparams.total_particles);

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
        loadCompoundatomsKernel << <boxparams.n_compounds, MAX_COMPOUND_PARTICLES >> > (atoms_dev, step, positions_dev, compounds_dev, coloringMethod, static_cast<float>(boxparams.boxSize));
    }
    if (boxparams.n_solvents > 0) {
        loadSolventatomsKernel << < boxparams.n_solvents / 128 + 1, 128 >> > (positions_dev, boxparams.n_compounds, boxparams.n_solvents, atoms_dev, static_cast<float>(boxparams.boxSize));   // TODO: This nsol/128 is wrong, if its not a multiple of 128
    }
}

const std::vector<RenderAtom>& Rasterizer::render(const Float3* positions, const std::vector<Compound>& compounds, const BoxParams& boxparams, int64_t step, Float3 camera_normal, ColoringMethod coloringMethod) {

    if (!isInitialized) { initialize(boxparams); }

    LIMA_UTILS::genericErrorCheck("Error before renderer");
	getAllAtoms(positions, compounds, boxparams, step, coloringMethod);

    cudaMemcpy(renderAtomsHost.data(), atoms_dev, sizeof(RenderAtom) * boxparams.total_particles, cudaMemcpyDeviceToHost);

    std::sort(renderAtomsHost.begin(), renderAtomsHost.end(), [camera_normal](const RenderAtom& a, const RenderAtom& b) {
        float dotA = a.position.x * camera_normal.x + a.position.y * camera_normal.y + a.position.z * camera_normal.z;
        float dotB = b.position.x * camera_normal.x + b.position.y * camera_normal.y + b.position.z * camera_normal.z;

        // Compare dot products for sorting
        return !a.IsDisabled() && dotA > dotB;
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
    case 'S':
        return ATOM_TYPE::S;
    case 'l':
		return ATOM_TYPE::LIMA_CUSTOM;
    default:
        printf("Unknown atom type: %c\n", atom);
        return ATOM_TYPE::NONE;
    }
}

__device__ float4 getColor(ATOM_TYPE atom_type) {
    switch (atom_type)
    {
    case ATOM_TYPE::SOL:
        return float4{0x03 / 255.0f, 0xa9 / 255.0f, 0xf4 / 255.0f, 0.0f };
    case ATOM_TYPE::H:
        return float4{0xF1 / 255.0f, 0xF1 / 255.0f, 0xF1 / 255.0f, 0.0f };
    case ATOM_TYPE::O:
        return float4{0xE0 / 255.0f, 0x20 / 255.0f, 0x20 / 255.0f, 0.0f };
    case ATOM_TYPE::C:
        return float4{0x30 / 255.0f, 0x10 / 255.0f, 0x90 / 255.0f, 0.0f };
    case ATOM_TYPE::P:
        return float4{0xFC / 255.0f, 0xF7 / 255.0f, 0x5E / 255.0f, 0.0f };
    case ATOM_TYPE::N:
        return float4{0x2E / 255.0f, 0x8B / 255.0f, 0x57 / 255.0f, 0.0f };
    case ATOM_TYPE::S:
        return float4{0xF4 / 255.0f, 0xC4 / 255.0f, 0x30 / 255.0f, 0.0f };
    case ATOM_TYPE::NONE:
        return float4{0xFF / 255.0f, 0x00 / 255.0f, 0xFF / 255.0f, 0.0f };
    default:
        return float4{0.0f, 0.0f, 0.0f, 0.0f };
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
    case ATOM_TYPE::S:
        return 0.189f;
    case ATOM_TYPE::LIMA_CUSTOM:
	    return 0.1f;
    case ATOM_TYPE::NONE:
        return 3.f;
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

    const int outputIndex = compounds[compound_id].absoluteIndexOfFirstParticle + threadIdx.x;

    const Compound* compound = &compounds[compound_id];
    
    if (local_id < compound->n_particles) {

        RenderAtom atom{};
        const Float3 positionNormalized = positions[global_id] / boxLenNM - 0.5f;// normalize from -0.5->0.5
        atom.position = float4{ positionNormalized.x, positionNormalized.y, positionNormalized.z, 0.f };        

        ATOM_TYPE atomType = RAS_getTypeFromAtomletter(compound->atomLetters[local_id]);

        atom.position.w = getRadius(atomType) / boxLenNM * 4.f;

        if (coloringMethod == ColoringMethod::Atomname)
            atom.color = getColor(atomType);
        else if (coloringMethod == ColoringMethod::Charge) {
            const float chargeNormalized = (compound->atom_charges[local_id] + elementaryChargeToCoulombPerMole)  / (elementaryChargeToCoulombPerMole*2.f);
            atom.color = float4{ chargeNormalized, 0.f, (1.f - chargeNormalized), 0.f };
        }

        atoms[outputIndex] = atom;
    }
}

__global__ void loadSolventatomsKernel(const Float3* positions, int n_compounds, int n_solvents, RenderAtom* atoms, float boxLenNM)
{
    const int solvent_index = blockIdx.x * blockDim.x + threadIdx.x;
    const int particle_index = n_compounds * MAX_COMPOUND_PARTICLES + solvent_index;

    if (solvent_index < n_solvents) {

		RenderAtom atom{};
        const Float3 positionNormalized = positions[particle_index] / boxLenNM - 0.5f;// normalize from -0.5->0.5
        atom.position = float4{ positionNormalized.x, positionNormalized.y, positionNormalized.z, 0.f };

        //atom.pos = positions[particle_index];
        //atom.pos = atom.pos / boxLenNM - 0.5f;    // normalize from -0.5->0.5

		//atom.atom_type = SOL;
        atom.position.w = getRadius(ATOM_TYPE::SOL) / boxLenNM * 4.f;
		// Debug
		//float velocity = (atom.pos - SolventCoord{ solventblock_prev->origo, solventblock_prev->rel_pos[threadIdx.x] }.getAbsolutePositionLM()).len();
  //      const float velocity = 1.f;
  //      const float point1nm = NANO_TO_LIMA * 0.1f;
		//const float color_scalar = velocity / point1nm * 255.f;
		//const uint8_t color_red = static_cast<uint8_t>(cuda::std::__clamp_to_integral<uint8_t, float>(color_scalar));
		//atom.color = Int3(color_red, 0, 255 - color_red);



        //atom.color = Int3(0, 0, 255);
        atom.color = float4{ 0,0,1, 0 };

        atoms[particle_index] = atom;
    }
}
