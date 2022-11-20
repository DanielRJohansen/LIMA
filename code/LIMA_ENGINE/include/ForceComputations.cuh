#pragma once



__device__ Float3 calcLJForce(Float3* pos0, Float3* pos1, float* data_ptr, float* potE, float sigma, float epsilon, int type1 = -1, int type2 = -1) {
	// Calculates LJ force on p0	(attractive to p1. Negative values = repulsion )//
	// input positions in cartesian coordinates [nm]
	// sigma [nm]
	// epsilon [J/mol]->[(kg*nm^2)/(ns^2*mol)]
	// Returns force in J/mol*M		?????????????!?!?//

	// Directly from book
	float dist_sq = (*pos1 - *pos0).lenSquared();
	float s = sigma * sigma / dist_sq;								// [nm^2]/[nm^2] -> unitless
	s = s * s * s;
	float force_scalar = 24.f * epsilon * s / dist_sq * (1.f - 2.f * s);	// Attractive. Negative, when repulsive		[(kg*nm^2)/(nm^2*ns^2*mol)] ->----------------------	[(kg)/(ns^2*mol)]	


	*potE += 4 * epsilon * s * (s - 1.f);

#ifdef LIMA_VERBOSE
	Float3 ddd = (*pos1 - *pos0) * force_scalar;
	if (ddd.x != ddd.x) {
		printf("\nError is here\n");
		printf("Scalar %f dist_sq %f\n", force_scalar, dist_sq);
		pos0->print('0');
		pos0->print('1');
	}
	//if (threadIdx.x == 32 && blockIdx.x == 28 && force < -600000.f && force > -700000) {
	//if (abs(force) > 20e+8 && type1 == 200 && type2 == 1) {
	if (abs(force_scalar) > 20e+9 || *potE != *potE) {
		//printf("\nDist %f   D2 %f    force %.1f [MN]     sigma: %.3f \tepsilon %.0f  t1 %d t2 %d\n", (float)sqrt(dist_sq), (float) dist_sq, (float)force_scalar *1e-6, sigma, epsilon, type1, type2);
		//pos0->print('1');
		//pos1->print('2');
		//printf("Block %d Thread %d\n", blockIdx.x, threadIdx.x);
		//printf("Thread %d Block %d self %f %f %f other %f %f %f\n", threadIdx.x, blockIdx.x, pos0->x, pos0->y, pos0->z, pos1->x, pos1->y, pos1->z);
		//printf("Force %f %f %f\n", force_unit_vector.x * force, force_unit_vector.y * force, force_unit_vector.z * force);
	}

	{
		//Float3 force_vec = force_unit_vector * force;
		/*if (force_vec.x != force_vec.x) {
			//printf("Force: %f\n", force);
			//force_unit_vector.print('u');
			printf("Thread %d Block %d self %f %f %f other %f %f %f\n", threadIdx.x, blockIdx.x, pos0->x, pos0->y, pos0->z, pos1->x, pos1->y, pos1->z);

		}*/
	}
#endif
	return (*pos1 - *pos0) * force_scalar;										// GN/mol [(kg*nm)/(ns^2*mol)]
}