#pragma once

#include "LimaTypes.cuh"

class Trajectory {
	std::vector<Float3> data;	

public:

	const int nFrames;
	const int nAtoms;
	const Float3 boxSize;// nm
	const float lambda; // [ps] deltaTime

	Trajectory(int nFrames, int nAtoms, Float3 boxSize, float dt) 
		: nFrames(nFrames), nAtoms(nAtoms), boxSize(boxSize), lambda(dt*NANO_TO_PICO) {
		data.resize(nFrames * nAtoms);
	}

	void Set(int frame, int atom, const Float3& pos) {
		if (frame >= nFrames)
			throw std::runtime_error("Frame index out of bounds");
		if (atom >= nAtoms)
			throw std::runtime_error("Atom index out of bounds");

		data[frame * nAtoms + atom] = pos;
	}
	
	// Purposefully not const, as the trr api requires non-const pointers
	Float3* GetFrame(int frame) {
		if (frame >= nFrames)
			throw std::runtime_error("Frame index out of bounds");

		return data.data() + frame * nAtoms;
	}
};