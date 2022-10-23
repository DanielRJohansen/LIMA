#include "Bodies.cuh"


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