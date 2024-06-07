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

struct LJParameter : public BondtypeBase<1> {
	static const int n_atoms = 1;


	LJParameter(const std::array<std::string, n_atoms>& typenames) : BondtypeBase(typenames) {}
	LJParameter(const std::array<std::string, n_atoms>& typenames, ForceField_NB::ParticleParameters parameters) : BondtypeBase(typenames), parameters(parameters) {}


	ForceField_NB::ParticleParameters parameters;
};

struct Singlebondtype : public BondtypeBase<2>{
	static const int n_atoms = 2;
	Singlebondtype(const std::array<std::string, n_atoms>& typenames, float b0 = 0.f, float kb = 0.f) : BondtypeBase(typenames), b0(b0), kb(kb) {
		sort();
	}
	Singlebondtype(const std::array<int, n_atoms>& ids, const std::array<std::string, n_atoms>& typenames)
		: BondtypeBase(ids, typenames) {
		sort();
	}

	float b0{};	// [nm]
	float kb{};	// [J/(mol*nm^2)]


	void flip() {
		swap(bonded_typenames[0], bonded_typenames[1]);
	}

	void sort() override {
		if (bonded_typenames[1] < bonded_typenames[0]) {
			flip();
		}			
	}
	
	SingleBond ToStandardBondRepresentation() const {
		return SingleBond{ {0,0}, b0 * NANO_TO_LIMA, kb / (NANO_TO_LIMA * NANO_TO_LIMA) };	// *puke*.. I really need a better solution, than having 3 representations of the same fucking data
	}
};




struct Anglebondtype : public BondtypeBase<3> {
	static const int n_atoms = 3;
	Anglebondtype(const std::array<std::string, n_atoms>& typenames, float t0=0.f, float kt=0.f)
		: BondtypeBase(typenames), theta0(t0), ktheta(kt) 
	{
		sort();
	}
	Anglebondtype(const std::array<int, n_atoms>& ids, const std::array<std::string, n_atoms>& typenames)
		: BondtypeBase(ids, typenames) {
		sort();
	}

	float theta0{};	// [rad]
	float ktheta{}; // [J/mole/rad^2]

	void flip() {
		swap(bonded_typenames[0], bonded_typenames[2]);
	}

	void sort() override {
		if (bonded_typenames[2] < bonded_typenames[0]) {
			flip();
		}
	}
	AngleBond ToStandardBondRepresentation() const {
		return AngleBond{ {0,0, 0}, theta0, ktheta};	// *puke*.. I really need a better solution, than having 3 representations of the same fucking data
	}
};

struct Dihedralbondtype : public BondtypeBase<4> {
	static const int n_atoms = 4;
	Dihedralbondtype(const std::array<std::string, n_atoms>& typenames, float phi0=0.f, float kphi=0.f, int n=0) 
		: BondtypeBase(typenames), phi0(phi0), kphi(kphi), n(n) 
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

	float phi0{}; // [rad]
	float kphi{}; // [J/mole]
	int n{};


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
	DihedralBond ToStandardBondRepresentation() const {
		return DihedralBond{ {0,0, 0, 0}, phi0, kphi, static_cast<float>(n)}; // Cast here????????	// *puke*.. I really need a better solution, than having 3 representations of the same fucking data
	}
};

struct Improperdihedralbondtype : public BondtypeBase<4> {
	static const int n_atoms = 4;
	// i j k l - https://manual.gromacs.org/current/reference-manual/functions/bonded-interactions.html
	
	Improperdihedralbondtype(const std::array<std::string, n_atoms>& typenames, float psi0=0.f, float kpsi=0.f)
		: BondtypeBase(typenames), psi0(psi0), kpsi(kpsi)
	{
		sort();
	}
	Improperdihedralbondtype(const std::array<int, n_atoms>& ids, const std::array<std::string, n_atoms>& typenames)
		: BondtypeBase(ids, typenames) {
		sort();
	}

	float psi0{};
	float kpsi{};

	void flip() {
		std::swap(bonded_typenames[0], bonded_typenames[3]);
		std::swap(bonded_typenames[1], bonded_typenames[2]);
	}

	//TODO: Check with Ali that this is okay?!
	void sort() override {
		// Improper dihedrals cannot be sorted, as they are asymmetric
	}
	ImproperDihedralBond ToStandardBondRepresentation() const {
		return ImproperDihedralBond{ {0,0, 0, 0}, psi0, kpsi};	// *puke*.. I really need a better solution, than having 3 representations of the same fucking data
	}
};


