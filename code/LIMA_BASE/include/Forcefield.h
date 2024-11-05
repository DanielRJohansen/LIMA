#pragma once

#include <string>

#include "Utilities.h"
#include "LimaTypes.cuh"
#include "MDFiles.h"
#include <span>



/////
// These are the full representation of what information a AtomType or Bondtype contains.
// Instances of such bonds are defined more compactly in bodies.cuh
/////

struct AtomType {
	std::string name{};
	int atNum{};
	ForceField_NB::ParticleParameters parameters{};
	float charge{}; // [kilo C/mol]
	float mass{}; // [kg/mol]
	char ptype{};
};

struct SinglebondType : public SingleBond {
	std::array<std::string, 2> bonded_typenames; // i,j
	int func{};
	Parameters params;
};

struct AnglebondType : public AngleUreyBradleyBond {
	std::array<std::string, 3> bonded_typenames; // i,j,k
	int func{};
	Parameters params;
};

struct DihedralbondType : public DihedralBond {
	std::array<std::string, 4> bonded_typenames; // ijkl
	int func{};
	Parameters params;
};

struct ImproperDihedralbondType : public ImproperDihedralBond {
	std::array<std::string, 4> bonded_typenames; // i j k l - https://manual.gromacs.org/current/reference-manual/functions/bonded-interactions.html
	int func{};
	Parameters params;
};


///-----------------------------------------------------------------------------------///


class AtomtypeDatabase;

template <typename GenericBondType>
class ParameterDatabase;

class LIMAForcefield {
public:
	LIMAForcefield();
	LIMAForcefield(const GenericItpFile& file);
	LIMAForcefield(const LIMAForcefield&) = delete;
	~LIMAForcefield();

	int GetActiveLjParameterIndex(const std::string& query);
	ForceField_NB GetActiveLjParameters();

	int GetActiveTinymoltypeIndex(const std::string& query);
	ForcefieldTinymol GetTinymolTypes();
	

	template<typename GenericBond>
	const std::vector<typename GenericBond::Parameters>& GetBondParameters(const auto& query);

private:
	std::unique_ptr<AtomtypeDatabase> ljParameters;
	std::unique_ptr<AtomtypeDatabase> tinymolTypes;

	std::unique_ptr<ParameterDatabase<SinglebondType>> singlebondParameters;
	std::unique_ptr<ParameterDatabase<AnglebondType>> anglebondParameters;
	std::unique_ptr<ParameterDatabase<DihedralbondType>> dihedralbondParameters;
	std::unique_ptr<ParameterDatabase<ImproperDihedralbondType>> improperdihedralbondParameters;

	template<typename GenericBond>
	const std::vector<typename GenericBond::Parameters>& _GetBondParameters(const auto& query);

	void LoadFileIntoForcefield(const GenericItpFile& file);
};