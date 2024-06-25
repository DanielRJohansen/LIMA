
#include "Rasterizer.cuh"
#include "Utilities.h"

__global__ void loadCompoundatomsKernel(RenderAtom* atoms, const int step, const Float3* positions, const Compound* compounds, ColoringMethod coloringMethod, float boxLenNM);
__global__ void loadSolventatomsKernel(const Float3* positions, int n_compounds, int n_solvents, RenderAtom* atoms, float boxLenNM, int nCompoundParticles);

const bool drawSolvent = false;
const bool drawHydrogens = false;

void Rasterizer::initialize(const BoxParams& boxparams, const std::vector<Compound>& compounds) {
    cudaMallocManaged(&positions_dev, sizeof(Float3) * boxparams.total_particles_upperbound);
    cudaMallocManaged(&compounds_dev, sizeof(Compound) * boxparams.n_compounds);

    cudaMemcpy(compounds_dev, compounds.data(), sizeof(Compound) * boxparams.n_compounds, cudaMemcpyHostToDevice);

	isInitialized = true;
}

Rasterizer::~Rasterizer() {
    if (isInitialized) {
        cudaFree(positions_dev);
        cudaFree(compounds_dev);
    }
}

/// <summary>	/// Returns a pointer to a list of atoms on the device	/// </summary>
void Rasterizer::getAllAtoms(const Float3* positions, const std::vector<Compound>& compounds, const BoxParams& boxparams, int64_t step, ColoringMethod coloringMethod, RenderAtom* renderAtoms) {
    cudaMemcpy(positions_dev, positions, sizeof(Float3) * boxparams.total_particles_upperbound, cudaMemcpyHostToDevice);
    LIMA_UTILS::genericErrorCheck("Error before 1");

    if (boxparams.n_compounds > 0) {
        loadCompoundatomsKernel << <boxparams.n_compounds, MAX_COMPOUND_PARTICLES >> > (renderAtoms, step, positions_dev, compounds_dev, coloringMethod, static_cast<float>(boxparams.boxSize));
    }    LIMA_UTILS::genericErrorCheck("Error before 2");

    if (boxparams.n_solvents > 0) {
        loadSolventatomsKernel << < boxparams.n_solvents / 128 + 1, 128 >> > (positions_dev, boxparams.n_compounds, boxparams.n_solvents, renderAtoms, static_cast<float>(boxparams.boxSize), boxparams.total_compound_particles);   // TODO: This nsol/128 is wrong, if its not a multiple of 128
    }    LIMA_UTILS::genericErrorCheck("Error before 3");

}

void Rasterizer::render(const Float3* positions, const std::vector<Compound>& compounds, const BoxParams& boxparams, 
    int64_t step, Float3 camera_normal, ColoringMethod coloringMethod, RenderAtom* renderAtoms) {

    if (!isInitialized) { initialize(boxparams, compounds); }

    LIMA_UTILS::genericErrorCheck("Error before renderer");
	getAllAtoms(positions, compounds, boxparams, step, coloringMethod, renderAtoms);

    LIMA_UTILS::genericErrorCheck("Error after renderer\n");
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
        return float4{0x03 / 255.0f, 0xa9 / 255.0f, 0xf4 / 255.0f, 1.0f };
    case ATOM_TYPE::H:
        return float4{0xF1 / 255.0f, 0xF1 / 255.0f, 0xF1 / 255.0f, 1.0f };
    case ATOM_TYPE::O:
        return float4{0xE0 / 255.0f, 0x20 / 255.0f, 0x20 / 255.0f, 1.0f };
    case ATOM_TYPE::C:
        return float4{0x30 / 255.0f, 0x10 / 255.0f, 0x90 / 255.0f, 1.0f };
    case ATOM_TYPE::P:
        return float4{0xFC / 255.0f, 0xF7 / 255.0f, 0x5E / 255.0f, 1.0f };
    case ATOM_TYPE::N:
        return float4{0x2E / 255.0f, 0x8B / 255.0f, 0x57 / 255.0f, 1.0f };
    case ATOM_TYPE::S:
        return float4{0xF4 / 255.0f, 0xC4 / 255.0f, 0x30 / 255.0f, 1.0f };
    case ATOM_TYPE::NONE:
        return float4{0xFF / 255.0f, 0x00 / 255.0f, 0xFF / 255.0f, 1.0f };
    default:
        return float4{0.0f, 0.0f, 0.0f, 1.0f };
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
            atom.color = float4{ chargeNormalized, 0.f, (1.f - chargeNormalized), 1.f };
        }

        if (atomType == ATOM_TYPE::H && !drawHydrogens) {
			atom.position.w = 0.f;
		}

        atoms[outputIndex] = atom;
    }
}

__global__ void loadSolventatomsKernel(const Float3* positions, int n_compounds, int n_solvents, RenderAtom* atoms, float boxLenNM, int nCompoundParticles)
{
    const int solvent_index = blockIdx.x * blockDim.x + threadIdx.x;
    const int particle_index = n_compounds * MAX_COMPOUND_PARTICLES + solvent_index;
    const int outputIndex = nCompoundParticles + solvent_index;

    if (solvent_index < n_solvents) {

		RenderAtom atom{};
        const Float3 positionNormalized = positions[particle_index] / boxLenNM - 0.5f;// normalize from -0.5->0.5
        atom.position = float4{ positionNormalized.x, positionNormalized.y, positionNormalized.z, 0.f };
        atom.position.w = getRadius(ATOM_TYPE::SOL) / boxLenNM * 4.f;

        atom.color = float4{ 0,0,1, drawSolvent };

        atoms[outputIndex] = atom;
    }
}
