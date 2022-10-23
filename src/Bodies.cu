#include "Bodies.cuh"




__host__ void Compound::init() {
	center_of_mass = calcCOM();
	//printf("")
	//radius = singlebonds[0].reference_dist * n_particles * 0.5f;
	//center_of_mass.print('C');
	//printf("Radius %f\n", radius);
}

__host__ void Compound::initBondedLUT()
{
	for (int i = 0; i < MAX_COMPOUND_PARTICLES; i++) {
		for (int ii = 0; ii < MAX_COMPOUND_PARTICLES; ii++) {
			bondedparticles_lookup[i][ii] = 0;
		}
	}
	
}

__host__ Float3 Compound::calcCOM() {
	Float3 com;
	for (int i = 0; i < n_particles; i++) {
		com += (prev_positions[i] * (1.f / (float)n_particles));
		//com += (particles[i].pos_tsub1 * (1.f / (float) n_particles));
	}

	return com;
}

__host__ void Compound::addParticle(int atomtype_id, Float3 pos) {
	if (n_particles == MAX_COMPOUND_PARTICLES) {
		printf("ERROR: Cannot add particle to compound!\n");
		exit(1);
	}

	atom_types[n_particles] = atomtype_id;
	prev_positions[n_particles] = pos;
	n_particles++;
}

__host__ void Compound::addParticle(int atomtype_id, Float3 pos, int atomtype_color_id, int global_id) {
	if (n_particles == MAX_COMPOUND_PARTICLES) {
		printf("ERROR: Cannot add particle to compound!\n");
		exit(1);
	}

	atom_types[n_particles] = atomtype_id;
	prev_positions[n_particles] = pos;
	atom_color_types[n_particles] = atomtype_color_id;
#ifdef LIMA_DEBUGMODE
	particle_global_ids[n_particles] = global_id;
	//printf("%d global id\n", global_id);
#endif
	n_particles++;
}

__host__ void Compound::calcParticleSphere() {
	Float3 com = calcCOM();// calcCOM(compound);

	//float furthest = LONG_MIN;
	float furthest = FLT_MIN;
	//float closest = LONG_MAX;
	float closest = FLT_MAX;
	int closest_index = 0;

	for (int i = 0; i < n_particles; i++) {
		float dist = (prev_positions[i] - com).len();
		closest_index = dist < closest ? i : closest_index;
		closest = std::min(closest, dist);
		furthest = std::max(furthest, dist);
	}

	key_particle_index = closest_index;
	confining_particle_sphere = furthest;
}

//__device__ void Compound::loadMeta(Compound* compound) {
//	n_particles = compound->n_particles;
//	n_singlebonds = compound->n_singlebonds;
//	n_anglebonds = compound->n_anglebonds;
//	n_dihedrals = compound->n_dihedrals;
//}
//__device__ void Compound::loadData(Compound* compound) {
//	if (threadIdx.x < n_particles) {
//		prev_positions[threadIdx.x] = compound->prev_positions[threadIdx.x];
//		atom_types[threadIdx.x] = compound->atom_types[threadIdx.x];
//		lj_ignore_list[threadIdx.x] = compound->lj_ignore_list[threadIdx.x];
//		forces[threadIdx.x] = compound->forces[threadIdx.x];
//		compound->forces[threadIdx.x] = Float3(0.f);
//		//#ifdef LIMA_DEBUGMODE
//		particle_global_ids[threadIdx.x] = compound->particle_global_ids[threadIdx.x];
//		//#endif
//	}
//	else {
//		prev_positions[threadIdx.x] = Float3(-1.f);
//		atom_types[threadIdx.x] = 0;
//		//lj_ignore_list[threadIdx.x] = compound->lj_ignore_list[threadIdx.x];
//		forces[threadIdx.x] = Float3(0.f);
//		particle_global_ids[threadIdx.x] = 0;
//	}
//	for (int i = 0; (i * blockDim.x) < n_singlebonds; i++) {
//		int index = i * blockDim.x + threadIdx.x;
//		if (index < n_singlebonds)
//			singlebonds[index] = compound->singlebonds[index];
//	}
//	for (int i = 0; (i * blockDim.x) < n_anglebonds; i++) {
//		int index = i * blockDim.x + threadIdx.x;
//		if (index < n_anglebonds)
//			anglebonds[index] = compound->anglebonds[index];
//	}
//	for (int i = 0; (i * blockDim.x) < n_dihedrals; i++) {
//		int index = i * blockDim.x + threadIdx.x;
//		if (index < n_dihedrals)
//			dihedrals[index] = compound->dihedrals[index];
//	}
//}




/*
Molecule::Molecule() {	// Always returns a h2o molecule rn
	n_atoms = 3;
	atoms = new Atom[3];
	uint8_t white[3] = { 250, 250, 250 };
	uint8_t red[3] = { 250, 0, 0 };
	uint8_t blue[3] = { 0,0,250 };
	// Using Van der Waals radii.... https://en.m.wikipedia.org/wiki/Van_der_Waals_radius
	atoms[0] = Atom(Float3(0, 0, 0), 0.152, 15.999, red);
	atoms[1] = Atom(Float3(0, 0, -95.7 / 1000.f * 1.5), 0.110, 1.008, white);
	atoms[2] = Atom(Float3(0, 92.65 / 1000, 23.96 / 1000.f), 0.110, 1.008, white);// was 0.53



	

	// Calc com
	CoM = Float3(0, 0, 0);
	double accumulated_mass = 0;
	for (int i = 0; i < 3; i++) {
		CoM = CoM + (atoms[i].pos * atoms[i].mass);
		accumulated_mass += atoms[i].mass;
	}
	CoM = CoM * (1 / accumulated_mass);

	for (int i = 0; i < 3; i++) {
		atoms[i].pos = atoms[i].pos - CoM;
		//printf("Atom %d pos: %f %f %f\n", i, atoms[i].pos.x, atoms[i].pos.y, atoms[i].pos.z);
	}
	

	// Rotate molecule so first atom is directly on top of the CoM
	Float3 pitch_yaw_roll(
		asin(atoms[0].pos.y / atoms[0].pos.len()),
		asin(atoms[0].pos.x / atoms[0].pos.len()),
		0
	);
	Float3 molecule_normal = Float3(0, 0, 1)._rotateAroundOrigin(pitch_yaw_roll * -1);
	for (int i = 0; i < n_atoms; i++) {
		atoms[i].pos = atoms[i].pos.rotateAroundVector(pitch_yaw_roll, molecule_normal);
		//printf("Atom %d pos: %f %f %f\n",i,  atoms[i].pos.x, atoms[i].pos.y, atoms[i].pos.z);
	}	
	//printf("molecule normal: %f %f %f\n", molecule_normal.x, molecule_normal.y, molecule_normal.z);
	//printf("pyr: %f %f %f\n", pitch_yaw_roll.x, pitch_yaw_roll.y, pitch_yaw_roll.z);
	
	// Noncritical rotation to make it look nice
	for (int i = 0; i < n_atoms; i++) {
		atoms[i].pos = atoms[i].pos.rotateAroundVector(Float3(0, PI/2.f, 0), Float3(0,1,0));
	}

}

	*/

CompoundBridgeBundleCompact::CompoundBridgeBundleCompact(CompoundBridgeBundle* bundle) {
	n_bridges = bundle->n_bridges;
	//printf("Transferring %d %d bridges\n", n_bridges, bundle->n_bridges);
	for (int i = 0; i < n_bridges; i++) {
		compound_bridges[i] = CompoundBridgeCompact(&bundle->compound_bridges[i]);
		//printf("bridge %d has %d particles\n\n", i, compound_bridges[i].n_particles);
	}
}

Molecule::Molecule() {
	compounds = new Compound[MAX_COMPOUNDS];
	//compound_bridge_bundle = new CompoundBridgeBundleCompact;
}

Float3 Molecule::calcCOM() {
	Float3 com(0.f);
	for (int i = 0; i < n_compounds; i++) {
		com += (compounds[i].calcCOM() * (1.f / (float)n_compounds));
	}
	return com;
}