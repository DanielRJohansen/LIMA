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


class LjParameterDatabase;

template <typename GenericBondType>
class ParameterDatabase;

class LIMAForcefield {
public:
	LIMAForcefield();
	LIMAForcefield(const fs::path& path, std::shared_ptr<std::vector<AtomType>> activeLJParamtypes);
	LIMAForcefield(const LIMAForcefield&) = delete;
	~LIMAForcefield();


	int GetActiveLjParameterIndex(const std::string& query);


	template<typename GenericBond>
	const std::vector<typename GenericBond::Parameters>& GetBondParameters(const auto& query);

	const fs::path path;

private:
	std::unique_ptr<LjParameterDatabase> ljParameters;

	std::unique_ptr<ParameterDatabase<SinglebondType>> singlebondParameters;
	std::unique_ptr<ParameterDatabase<AnglebondType>> anglebondParameters;
	std::unique_ptr<ParameterDatabase<DihedralbondType>> dihedralbondParameters;
	std::unique_ptr<ParameterDatabase<ImproperDihedralbondType>> improperdihedralbondParameters;

	void LoadFileIntoForcefield(const fs::path& path);
};

class ForcefieldManager {
	std::shared_ptr<std::vector<AtomType>> activeLJParamtypes;

	std::vector<std::unique_ptr<LIMAForcefield>> forcefields;

	LIMAForcefield& GetForcefield(const fs::path& forcefieldName);	

	const fs::path internalForcefieldsDir = FileUtils::GetLimaDir() / "resources/forcefields";

	const fs::path limaTestForcefield = internalForcefieldsDir / "lima_custom_forcefield.itp";
	const fs::path defaultForcefield = internalForcefieldsDir / "charmm27.ff/forcefield.itp";

public:

	ForcefieldManager();
	~ForcefieldManager();

	int GetActiveLjParameterIndex(const std::vector<fs::path>& forcefieldName, const std::string& query);
	ForceField_NB GetActiveLjParameters();

	template<typename GenericBond>
	const std::vector<typename GenericBond::Parameters>& GetBondParameters(
		const std::vector<fs::path>& forcefieldNames, const auto& query);
};
