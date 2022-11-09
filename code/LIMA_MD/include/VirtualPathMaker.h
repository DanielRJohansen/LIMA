#pragma once

#include "LimaTypes.cuh"
#include <math.h>

class VirtualPathMaker
{
public:
	Float3* makeVirtualPath(Float3* particle_positions, Float3* forces, int n_steps);




private:
	Float3** generateAllPaths(Float3* particle_positions, Float3* forces, int n_steps);
	int findIndexOfShortestStep(Float3 pos_prev, Float3* possible_positions);
	int findIndexOfShortestPath(Float3** all_paths, int n_steps);


	Float3* generateAllPositions(Float3 particle_position, Float3 force);

	double getAtttractivePosition(double force);	// returns position relative to particle_position!!
	double getRepulsivePosition(double force);	// returns position relative to particle_position!!



	double binarySearch(double lower, double middle, double upper, double force);
	double calcForce(double distance);
	double lengthOfPath(Float3* path, int n_steps);


	const double sigma = 0.3923f;	//nm, basicllay break 0 point
	const double epsilon = 0.5986 * 1'000.f; // J/mol
	double threshold = 0.1; // N/mol
};

