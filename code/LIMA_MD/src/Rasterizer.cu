#include "LIMA_MD/src/Rasterizer.cuh"
#include "LIMA_ENGINE/include/EngineUtils.cuh"
#include "LIMA_BASE/include/Utilities.h"


#include <algorithm>
#include <cuda/std/cmath>



std::vector<RenderBall> Rasterizer::render(Simulation* simulation) {
	solvent_offset = simulation->n_compounds * MAX_COMPOUND_PARTICLES;
	n_threadblocks = (int)ceil((float)simulation->total_particles_upperbound / (float)RAS_THREADS_PER_BLOCK);

    LIMA_UTILS::genericErrorCheck("Error before renderer");
	RenderAtom* atoms_dev = getAllAtoms(simulation);

    std::vector<RenderBall> balls_host = processAtoms(atoms_dev, simulation);
    cudaFree(atoms_dev);

    std::sort(balls_host.begin(), balls_host.end(), [](const RenderBall& a, const RenderBall& b) {return !a.disable && a.pos.y < b.pos.y; });    

    LIMA_UTILS::genericErrorCheck("Error after renderer");

	return balls_host;
}












__global__ void loadCompoundatomsKernel(Box* box, RenderAtom* atoms, const int step);
__global__ void loadSolventatomsKernel(Box* box, RenderAtom* atoms, int offset, const int step);
__global__ void processAtomsKernel(RenderAtom* atoms, RenderBall* balls);

__global__ void kernelA(Box* box) {
    // Do nothing
}
__global__ void kernelB(Box* box) {
    // Do nothing
}

RenderAtom* Rasterizer::getAllAtoms(Simulation* simulation) {
	RenderAtom* atoms;
	cudaMallocManaged(&atoms, sizeof(RenderAtom) * simulation->total_particles_upperbound);


    //kernelA << <10, 10>> > (simulation->sim_dev->box);
    //cudaDeviceSynchronize();
    //kernelB << <10, 10>> > (simulation->sim_dev->box);


    Box* boxptr = simulation->sim_dev->box; // Use intermediate ptr to avoid prepascal limitation of concurrent managed data access
    if (simulation->n_compounds > 0) {
        loadCompoundatomsKernel << <simulation->n_compounds, MAX_COMPOUND_PARTICLES >> > (boxptr, atoms, simulation->simparams_host.step);
    }
	if (simulation->n_solvents > 0) {
		loadSolventatomsKernel << < SolventBlockGrid::blocks_total, MAX_SOLVENTS_IN_BLOCK >> > (boxptr, atoms, solvent_offset, simulation->simparams_host.step);
	}

	return atoms;
}


std::vector<RenderBall> Rasterizer::processAtoms(RenderAtom* atoms, Simulation* simulation) {
    RenderBall* balls_device;
    cudaMalloc(&balls_device, sizeof(RenderBall) * simulation->total_particles_upperbound);
    processAtomsKernel <<< n_threadblocks, RAS_THREADS_PER_BLOCK >>> (atoms, balls_device);
    LIMA_UTILS::genericErrorCheck("Error during rendering");

    std::vector<RenderBall> balls_host(simulation->total_particles_upperbound);    
    cudaMemcpy(balls_host.data(), balls_device, sizeof(RenderBall) * simulation->total_particles_upperbound, cudaMemcpyDeviceToHost);

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

__device__ float getRadius(ATOM_TYPE atom_type) {
    switch (atom_type)
    {

    case ATOM_TYPE::H:
        return 0.04;
    case ATOM_TYPE::C:
        return 0.1;
    case ATOM_TYPE::N:
        return 0.1;
    case ATOM_TYPE::O:
        return 0.12;
    case ATOM_TYPE::SOL:
        return 0.05;
    case ATOM_TYPE::P:
        return 0.15;
    case ATOM_TYPE::NONE:
        return 1;
    default:
        return 1;
    }
}









__global__ void loadCompoundatomsKernel(Box* box, RenderAtom* atoms, const int step) {

    int local_id = threadIdx.x;
    int compound_id = blockIdx.x;
    int global_id = threadIdx.x + blockIdx.x * blockDim.x;

    Compound* compound = &box->compounds[compound_id];
    
    if (local_id < compound->n_particles) {
        auto coordarray_ptr = CoordArrayQueueHelpers::getCoordarrayRef(box->coordarray_circular_queue, step, compound_id);

        RenderAtom atom{};
        //atom.pos = LIMAPOSITIONSYSTEM::getGlobalPosition(*coordarray_ptr);
        atom.pos = LIMAPOSITIONSYSTEM::getAbsolutePositionNM(coordarray_ptr->origo, coordarray_ptr->rel_positions[local_id]);


        //atoms[global_id].pos.print('A');
        atom.mass = SOLVENT_MASS;                                                         // TEMP
        //atoms[global_id].atom_type = RAS_getTypeFromIndex(box->compounds[compound_id].atom_types[local_id]);
        atom.atom_type = RAS_getTypeFromIndex(compound->atom_color_types[local_id]);

        atom.color = getColor(atom.atom_type);

        atoms[global_id] = atom;
    }
    else {
        atoms[global_id].atom_type = ATOM_TYPE::NONE;
    }
}

__global__ void loadSolventatomsKernel(Box* box, RenderAtom* atoms, int offset, const int step) {
    SolventBlock* solventblock = CoordArrayQueueHelpers::getSolventBlockPtr(box->solventblockgrid_circular_queue, step, blockIdx.x);
    SolventBlock* solventblock_prev = CoordArrayQueueHelpers::getSolventBlockPtr(box->solventblockgrid_circular_queue, step == 0 ? 0 : step-1, blockIdx.x);

    if (threadIdx.x < solventblock->n_solvents) {
        //const SolventCoord coord{solventblock->origo, solventblock->rel_pos[threadIdx.x] };

		RenderAtom atom{};
		//atom.pos = coord.getAbsolutePositionLM();
        atom.pos = LIMAPOSITIONSYSTEM::getAbsolutePositionNM(solventblock->origo, solventblock->rel_pos[threadIdx.x]);
        EngineUtils::applyPBCNM(&atom.pos);   // TMP, dunno if i wanna do this.
		atom.mass = SOLVENT_MASS;
		atom.atom_type = SOL;

		// Debug
		//float velocity = (atom.pos - SolventCoord{ solventblock_prev->origo, solventblock_prev->rel_pos[threadIdx.x] }.getAbsolutePositionLM()).len();
        float velocity = 1.f;
        float point1nm = NANO_TO_LIMA * 0.1f;
		float color_scalar = velocity / point1nm * 255.f;
		uint8_t color_red = static_cast<uint8_t>(cuda::std::__clamp_to_integral<uint8_t, float>(color_scalar));
		atom.color = Int3(color_red, 0, 255 - color_red);
		//printf("vel %f, %d\n", velocity, color_red);

        // This part is for various debugging purposes
        //int query_id = 0;
        //if (solvent_id == query_id) {
        //    atoms[solvent_id + offset].atom_type = P;
        //}
        //const auto& nlist = box->solvent_neighborlists[solvent_id];
        //for (int i = 0; i < nlist.n_solvent_neighbors; i++) {
        //    if (nlist.neighborsolvent_ids[i] == query_id) {
        //        atoms[solvent_id + offset].atom_type = O;
        //    }
        //}
        int global_id = threadIdx.x + blockIdx.x * blockDim.x;
        atoms[global_id] = atom;
        //if (solvent_id != 440)
        //    atoms[solvent_id + offset].atom_type = NONE;
    }
}

__global__ void processAtomsKernel(RenderAtom* atoms, RenderBall* balls) { 
    const int index = threadIdx.x + blockIdx.x * RAS_THREADS_PER_BLOCK;
    
    RenderAtom atom = atoms[index];

    
    //atom.color = getColor(atom.atom_type);

    atom.radius = (getRadius(atom.atom_type)) / (1.f+atom.pos.y * 0.00000000001f);       // [nm]

    // Convert units to normalized units for OpenGL
    atom.radius = 0.25f * atom.radius;            // Yeah, i'm just eyeballing this..

    for (int dim = 0; dim < 3; dim++) {
        *atom.pos.placeAt(dim) = (atom.pos.at(dim) / BOX_LEN_NM - 0.5f) *1.8f;
    }
    
    const RenderBall ball(atom.pos, atom.radius, atom.color, atom.atom_type);
    balls[index] = ball;
}
