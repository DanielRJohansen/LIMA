#pragma once

#include "LimaTypes.cuh"
#include "MoleculeHull.cuh"

#include <functional>
#include <optional>

class ConvexHullEngine {
public:

	void MoveMoleculesUntillNoOverlap(MoleculeHullCollection& mhCol, Float3 boxSize, 
		std::function<void()> callEachIteration = []() {});



};