#pragma once

#include "Bodies.cuh"
#include "Utilities.h"

#include <algorithm>
#include <array>
#include <format>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

///////////////////////////////// READ HERE FIRST /////////////////////////////////
// ffbonded.itp and ffnonbonden.itp has different atomtypes. 
// The purpose if this bit of code is to compress all forcefield parameters so they can fit on a GPU more easily
//


//template <typename bond>
//void SortBondedtypeNames(std::array<std::string, bond::n_atoms>& bondedTypeNames) {
//	if constexpr (std::is_same<bond, SingleBond>::value) {
//		if (bondedTypeNames[1] < bondedTypeNames[0]) {
//			swap(bondedTypeNames[0], bondedTypeNames[1]);
//		}
//	}
//	else if constexpr (std::is_same<bond, AngleBond>::value) {
//		if (bondedTypeNames[2] < bondedTypeNames[0]) {
//			swap(bondedTypeNames[0], bondedTypeNames[2]);
//		}
//	}
//	else if constexpr (std::is_same<bond, DihedralBond>::value) {
//		// If out is wrong, flip
//		if (bondedTypeNames[3] < bondedTypeNames[0]) {
//			std::swap(bondedTypeNames[0], bondedTypeNames[3]);
//			std::swap(bondedTypeNames[1], bondedTypeNames[2]);
//		}
//		// If outer is identical, but inner is wrong, flip
//		else if (
//			bondedTypeNames[0] == bondedTypeNames[3] &&
//			bondedTypeNames[2] < bondedTypeNames[1])
//		{
//			std::swap(bondedTypeNames[0], bondedTypeNames[3]);
//			std::swap(bondedTypeNames[1], bondedTypeNames[2]);
//		}
//	}
//	else if constexpr (std::is_same<bond, ImproperDihedralBond>::value) {
//		if (bondedTypeNames[3] < bondedTypeNames[0]) {
//			std::swap(bondedTypeNames[0], bondedTypeNames[3]);
//		}
//	}
//	else {
//		throw std::runtime_error("Illegal bond type");
//	}
//}




struct AtomType {
	std::string name{};
	int atNum{};
	ForceField_NB::ParticleParameters parameters{};
	float charge{}; // [kilo C/mol]
	char ptype{};
};

struct SinglebondType {
	static const int nAtoms = 2;
	std::array<std::string, 2> bonded_typenames; // i,j
	int func{};
	SingleBond::Parameters params;
};

struct AnglebondType {
	static const int nAtoms = 3;
	std::array<std::string, 3> bonded_typenames; // i,j,k
	int func{};
	AngleBond::Parameters params;
	float ub0;	// No idea what this is
	float cub;	// No idea what this is
};

struct DihedralbondType {
	static const int nAtoms = 4;
	std::array<std::string, 4> bonded_typenames; // ijkl
	int func{};
	DihedralBond::Parameters params;
};

struct ImproperDihedralbondType {
	static const int nAtoms = 4;
	std::array<std::string, 4> bonded_typenames; // i j k l - https://manual.gromacs.org/current/reference-manual/functions/bonded-interactions.html
	int func{};
	ImproperDihedralBond::Parameters params;
};

template <typename bond>
void FlipBondedtypeNames(std::array<std::string, bond::nAtoms>& bondedTypeNames) {
	if constexpr (std::is_same<bond, SinglebondType>::value) {
		swap(bondedTypeNames[0], bondedTypeNames[1]);
	}
	else if constexpr (std::is_same<bond, AnglebondType>::value) {
		swap(bondedTypeNames[0], bondedTypeNames[2]);
	}
	else if constexpr (std::is_same<bond, DihedralbondType>::value) {

		std::swap(bondedTypeNames[0], bondedTypeNames[3]);
		std::swap(bondedTypeNames[1], bondedTypeNames[2]);

	}
	else if constexpr (std::is_same<bond, ImproperDihedralbondType>::value) {
		std::swap(bondedTypeNames[0], bondedTypeNames[3]);
	}
	else {
		throw std::runtime_error("Illegal bond type");
	}
}
