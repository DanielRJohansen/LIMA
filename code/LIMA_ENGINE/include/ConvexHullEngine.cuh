#pragma once

#include "LimaTypes.cuh"
#include "MoleculeHull.cuh"


class ConvexHullEngine {
public:

	void MoveMoleculesUntillNoOverlap(MoleculeHullCollection& mhCol, Float3 boxSize);


};