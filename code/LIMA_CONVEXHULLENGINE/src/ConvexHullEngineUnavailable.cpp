// Replaces ConvexHullEngine.cu in HIP builds, where it isn't ported yet
#include "ConvexHullEngine.cuh"

#include <stdexcept>

void ConvexHullEngine::MoveMoleculesUntillNoOverlap(MoleculeHullCollection&, Float3, std::function<void()>) {
	throw std::runtime_error("ConvexHullEngine is not available in HIP builds");
}
