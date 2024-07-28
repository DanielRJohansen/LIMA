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


template <typename bond>
void SortBondedtypeNames(std::array<std::string, bond::n_atoms>& bondedTypeNames) {
	if constexpr (std::is_same<bond, SingleBond>::value) {
		if (bondedTypeNames[1] < bondedTypeNames[0]) {
			swap(bondedTypeNames[0], bondedTypeNames[1]);
		}
	}
	else if constexpr (std::is_same<bond, AngleBond>::value) {
		if (bondedTypeNames[2] < bondedTypeNames[0]) {
			swap(bondedTypeNames[0], bondedTypeNames[2]);
		}
	}
	else if constexpr (std::is_same<bond, DihedralBond>::value) {
		// If out is wrong, flip
		if (bondedTypeNames[3] < bondedTypeNames[0]) {
			std::swap(bondedTypeNames[0], bondedTypeNames[3]);
			std::swap(bondedTypeNames[1], bondedTypeNames[2]);
		}
		// If outer is identical, but inner is wrong, flip
		else if (
			bondedTypeNames[0] == bondedTypeNames[3] &&
			bondedTypeNames[2] < bondedTypeNames[1])
		{
			std::swap(bondedTypeNames[0], bondedTypeNames[3]);
			std::swap(bondedTypeNames[1], bondedTypeNames[2]);
		}
	}
	else if constexpr (std::is_same<bond, ImproperDihedralBond>::value) {
		if (bondedTypeNames[3] < bondedTypeNames[0]) {
			std::swap(bondedTypeNames[0], bondedTypeNames[3]);
		}
	}
	else {
		throw std::runtime_error("Illegal bond type");
	}
}


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
	SingleBond::Parameters parameters;

	void flip() {
		swap(bonded_typenames[0], bonded_typenames[1]);
	}

	void sort() {
		if (bonded_typenames[1] < bonded_typenames[0]) {
			flip();
		}
	}
};

struct AnglebondType {
	static const int nAtoms = 3;
	std::array<std::string, 3> bonded_typenames; // i,j,k
	int func{};
	AngleBond::Parameters parameters;
	float ub0;	// No idea what this is
	float cub;	// No idea what this is

	void flip() {
		swap(bonded_typenames[0], bonded_typenames[2]);
	}
	void sort() {
		if (bonded_typenames[2] < bonded_typenames[0]) {
			flip();
		}
	}
};

struct DihedralbondType {
	static const int nAtoms = 4;
	std::array<std::string, 4> bonded_typenames; // ijkl
	int func{};
	DihedralBond::Parameters parameters;

	void flip() {
		std::swap(bonded_typenames[0], bonded_typenames[3]);
		std::swap(bonded_typenames[1], bonded_typenames[2]);
	}
	void sort() {
		// If out is wrong, flip
		if (bonded_typenames[3] < bonded_typenames[0]) {
			flip();
		}
		// If outer is identical, but inner is wrong, flip
		else if (
			bonded_typenames[0] == bonded_typenames[3] &&
			bonded_typenames[2] < bonded_typenames[1])
		{
			flip();
		}
	}
};

struct ImproperDihedralbondType {
	static const int nAtoms = 4;
	std::array<std::string, 4> bonded_typenames; // ijkl
	int func{};
	ImproperDihedralBond::Parameters parameters;

	void flip() {
		std::swap(bonded_typenames[0], bonded_typenames[3]);
		//std::swap(bonded_typenames[1], bonded_typenames[2]);
	}
	//TODO: Check with Ali that this is okay?!
	void sort() {
		if (bonded_typenames[3] < bonded_typenames[0]) {
			flip();
		}
		// Improper dihedrals cannot be sorted, as they are asymmetric
		//flip();	// TODO: I dont know if this needs to happen for all forcefields, but it does for CHARMM
	}
};


template <int n_atoms>	// n atoms in bond
struct BondtypeBase {
	BondtypeBase(const std::array<std::string, n_atoms>& typenames) : bonded_typenames(typenames) {
		for (int i = 0; i < n_atoms; i++) {
			global_ids[i] = -1;
		}
	}
	BondtypeBase(const std::array<int, n_atoms>& ids, const std::array<std::string, n_atoms>& typenames)
		: bonded_typenames(typenames), global_ids(ids) {
	}

	virtual void sort() {};
	virtual void flip() {};

	std::array<std::string, n_atoms> bonded_typenames;
	std::array<int, n_atoms> global_ids;
};

struct Singlebondtype : public BondtypeBase<2>{
	static const int n_atoms = 2;
	Singlebondtype(const std::array<std::string, n_atoms>& typenames, float b0 = 0.f, float kb = 0.f) 
		: BondtypeBase(typenames), params{ b0, kb } {
		sort();
	}
	Singlebondtype(const std::array<int, n_atoms>& ids, const std::array<std::string, n_atoms>& typenames)
		: BondtypeBase(ids, typenames) {
		sort();
	}

	SingleBond::Parameters params;


	void flip() {
		swap(bonded_typenames[0], bonded_typenames[1]);
	}

	void sort() override {
		if (bonded_typenames[1] < bonded_typenames[0]) {
			flip();
		}			
	}
	
	//SingleBond ToStandardBondRepresentation() const {
	//	return SingleBond{ {0,0}, b0 * NANO_TO_LIMA, kb / (NANO_TO_LIMA * NANO_TO_LIMA) };	// TODO: *puke*.. I really need a better solution, than having 3 representations of the same fucking data
	//}
};




struct Anglebondtype : public BondtypeBase<3> {
	static const int n_atoms = 3;
	Anglebondtype(const std::array<std::string, n_atoms>& typenames, float t0=0.f, float kt=0.f)
		: BondtypeBase(typenames), params{ t0, kt }
	{
		sort();
	}
	Anglebondtype(const std::array<int, n_atoms>& ids, const std::array<std::string, n_atoms>& typenames)
		: BondtypeBase(ids, typenames) {
		sort();
	}
	AngleBond::Parameters params;

	void flip() {
		swap(bonded_typenames[0], bonded_typenames[2]);
	}

	void sort() override {
		if (bonded_typenames[2] < bonded_typenames[0]) {
			flip();
		}
	}
	//AngleBond ToStandardBondRepresentation() const {
	//	return AngleBond{ {0,0, 0}, theta0, ktheta};	// *puke*.. I really need a better solution, than having 3 representations of the same fucking data
	//}
};

struct Dihedralbondtype : public BondtypeBase<4> {
	static const int n_atoms = 4;
	Dihedralbondtype(const std::array<std::string, n_atoms>& typenames, float phi0=0.f, float kphi=0.f, int n=0) 
		: BondtypeBase(typenames), params{phi0, kphi, n}
	{
		if (typenames[0] == "X" && typenames[3] == "X") {
			int a = 0;
		}
		sort();
	}
	Dihedralbondtype(const std::array<int, n_atoms>& ids, const std::array<std::string, n_atoms>& typenames)
		: BondtypeBase(ids, typenames) {
		sort();
	}
	DihedralBond::Parameters params;


	void flip() {
		std::swap(bonded_typenames[0], bonded_typenames[3]);
		std::swap(bonded_typenames[1], bonded_typenames[2]);
	}

	void sort() override {
		// If out is wrong, flip
		if (bonded_typenames[3] < bonded_typenames[0]) {
			flip();
		}
		// If outer is identical, but inner is wrong, flip
		else if (
			bonded_typenames[0] == bonded_typenames[3] && 
			bonded_typenames[2] < bonded_typenames[1]) 
		{
			flip();
		}
	}
	//DihedralBond ToStandardBondRepresentation() const {
	//	return DihedralBond{ {0,0, 0, 0}, phi0, kphi, static_cast<float>(n)}; // Cast here????????	// *puke*.. I really need a better solution, than having 3 representations of the same fucking data
	//}
};

struct Improperdihedralbondtype : public BondtypeBase<4> {
	static const int n_atoms = 4;
	// i j k l - https://manual.gromacs.org/current/reference-manual/functions/bonded-interactions.html
	
	Improperdihedralbondtype(const std::array<std::string, n_atoms>& typenames, float psi0=0.f, float kpsi=0.f)
		: BondtypeBase(typenames), params{ psi0, kpsi }
	{
		sort();
	}
	ImproperDihedralBond::Parameters params;

	void flip() {
		std::swap(bonded_typenames[0], bonded_typenames[3]);
		//std::swap(bonded_typenames[1], bonded_typenames[2]);
	}

	//TODO: Check with Ali that this is okay?!
	void sort() override {
		if (bonded_typenames[3] < bonded_typenames[0]) {
			flip();
		}
		// Improper dihedrals cannot be sorted, as they are asymmetric
		//flip();	// TODO: I dont know if this needs to happen for all forcefields, but it does for CHARMM
	}
	//ImproperDihedralBond ToStandardBondRepresentation() const {
	//	return ImproperDihedralBond{ {0,0, 0, 0}, psi0, kpsi};	// *puke*.. I really need a better solution, than having 3 representations of the same fucking data
	//}
};


