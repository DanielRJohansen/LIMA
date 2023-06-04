#include "LIMA_MD/src/Rasterizer.cuh"
#include "LIMA_ENGINE/include/EngineUtils.cuh"






void mergeSortAPI(RenderBall* balls, int n_balls);


RenderBall* Rasterizer::render(Simulation* simulation) {
	solvent_offset = simulation->n_compounds * MAX_COMPOUND_PARTICLES;
	n_threadblocks = (int)ceil((float)simulation->total_particles_upperbound / (float)RAS_THREADS_PER_BLOCK);


	RenderAtom* atoms = getAllAtoms(simulation);


	RenderBall* balls = processAtoms(atoms, simulation);
	mergeSortAPI(balls, simulation->total_particles_upperbound);

	return balls;
}











__global__ void loadCompoundatomsKernel(Box* box, RenderAtom* atoms, const int step);
__global__ void loadSolventatomsKernel(Box* box, RenderAtom* atoms, int offset, const int step);
__global__ void processAtomsKernel(RenderAtom* atoms, RenderBall* balls);

RenderAtom* Rasterizer::getAllAtoms(Simulation* simulation) {


	RenderAtom* atoms;
	cudaMalloc(&atoms, sizeof(RenderAtom) * simulation->total_particles_upperbound);
	//cudaMalloc(&atoms, sizeof(RenderAtom) * SolventBlockGrid::blocks_total * MAX_SOLVENTS_IN_BLOCK);


	//int solvent_blocks = (int)ceil((float)simulation->n_solvents / (float)THREADS_PER_LOADSOLVENTSATOMSKERNEL);


	Box* box = simulation->box;
	if (simulation->n_compounds > 0)
		loadCompoundatomsKernel << <simulation->n_compounds, MAX_COMPOUND_PARTICLES >> > (box, atoms, simulation->simparams_host.step);
	if (simulation->n_solvents > 0) {
		//loadSolventatomsKernel << < solvent_blocks, THREADS_PER_LOADSOLVENTSATOMSKERNEL >> > (simulation->box, atoms, solvent_offset);
		loadSolventatomsKernel << < SolventBlockGrid::blocks_total, MAX_SOLVENTS_IN_BLOCK >> > (simulation->box, atoms, solvent_offset, simulation->simparams_host.step);
	}

	cudaDeviceSynchronize();

	return atoms;
}

void Rasterizer::sortAtoms(RenderAtom* atoms, int dim) {

    cudaDeviceSynchronize();
}

RenderBall* Rasterizer::processAtoms(RenderAtom* atoms, Simulation* simulation) {
    RenderBall* balls_device;
    cudaMalloc(&balls_device, sizeof(RenderBall) * simulation->total_particles_upperbound);
    processAtomsKernel <<< n_threadblocks, RAS_THREADS_PER_BLOCK >>> (atoms, balls_device);
    cudaDeviceSynchronize();

    RenderBall* balls_host = new RenderBall[simulation->total_particles_upperbound];
    cudaMemcpy(balls_host, balls_device, sizeof(RenderBall) * simulation->total_particles_upperbound, cudaMemcpyDeviceToHost);

    cudaFree(atoms);
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









__global__ void loadCompoundatomsKernel(Box* box, RenderAtom* atoms, const int step) {                                                            // TODO: CAN ONLY HANDLE ONE COMPOUND!!!
    int local_id = threadIdx.x;
    int compound_id = blockIdx.x;
    int global_id = threadIdx.x + blockIdx.x * blockDim.x;

    
    if (local_id < box->compounds[compound_id].n_particles) {
        //atoms[global_id].pos = LIMAPOSITIONSYSTEM::getGlobalPosition(box->compound_coord_array[compound_id]);
        auto coordarray_ptr = CoordArrayQueueHelpers::getCoordarrayRef(box->coordarray_circular_queue, step, compound_id);

        RenderAtom atom{};
        //atom.pos = LIMAPOSITIONSYSTEM::getGlobalPosition(*coordarray_ptr);
        atom.pos = LIMAPOSITIONSYSTEM::getAbsolutePositionNM(coordarray_ptr->origo, coordarray_ptr->rel_positions[local_id]);


        //atoms[global_id].pos.print('A');
        atom.mass = SOLVENT_MASS;                                                         // TEMP
        //atoms[global_id].atom_type = RAS_getTypeFromIndex(box->compounds[compound_id].atom_types[local_id]);
        atom.atom_type = RAS_getTypeFromIndex(box->compounds[compound_id].atom_color_types[local_id]);

        atom.color = getColor(atom.atom_type);

        atoms[global_id] = atom;
    }
    else {
        atoms[global_id].atom_type = ATOM_TYPE::NONE;
    }
}
#include <algorithm> 
#include <cuda/std/cmath>

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


RenderBall* merge(const RenderBall* left, int n_left, const RenderBall* right, int n_right) {
    RenderBall* merged = new RenderBall[n_left + n_right];
    int l = 0;
    int r = 0;
    int index = 0;

    while (l < n_left && r < n_right) {
        if (left[l].pos.y < right[r].pos.y) {
            merged[index++] = left[l++];
        }
        else {
            merged[index++] = right[r++];
        }
    }
    while (l < n_left) { merged[index++] = left[l++]; }
    while (r < n_right) { merged[index++] = right[r++]; }


    return merged;
}

RenderBall* mergeSort(const RenderBall* atoms, const int l, const int r) {	// l and r are indexes of the two extremes
    int m = (l + r) / 2;

    RenderBall* left;
    RenderBall* right;

    //printf("step\n");
    if (r - l > 1) {
        //printf("l %d r %d\n", l, r);
        left = mergeSort(atoms, l, m - 1);
        right = mergeSort(atoms, m, r);
        RenderBall* merged = merge(left, m - l, right, r - m + 1);

        if (l == m - 1)		// Take care of special case, where only left side can be a single object instead of array!
            delete left;
        else
            delete[] left;
        delete[] right;

        return merged;
    }
    else if (r - l == 1) {
        return merge(&atoms[l], 1, &atoms[r], 1);
    }
    else {
        return new RenderBall(atoms[l]);
    }
}

void mergeSortAPI(RenderBall* balls, int n_balls) {					// Returns a mapping where mapping[0] is the closest id, mapping [1] seconds closest, so on


    RenderBall* sorted_balls = mergeSort(balls, 0, n_balls - 1);
    for (int i = 0; i < n_balls; i++) {
        balls[i] = sorted_balls[i];
    }

    int* mapping = new int[n_balls];
    //for (int i = 0; i < n_atoms; i++) mapping[i] = sorted_atoms[i].id;

    delete[] sorted_balls;
}