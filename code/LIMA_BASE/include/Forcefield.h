#pragma once

#include <string>

#include "Utilities.h"
#include "LimaTypes.cuh"
#include "MDFiles.h"
#include <span>
#include <optional>
#include <unordered_map>



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

struct PairbondType : public PairBond {
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
	std::vector<NonbondedInteractionParams> GetNonbondedInteractionParams() const;

	NBParams GetLjParameters(const std::string& query) const;// Experimental, breaks the ActiveLjParamsSystem
	std::optional<AtomType> GetAtomtype(const std::string& query) const;
//	ForcefieldTinymol GetTinymolTypes();
	

	template<typename GenericBond>
	const std::vector<typename GenericBond::Parameters>& GetBondParameters(const auto& query);
	// Periodic impropers (funct 4), which GROMACS matches only against [ dihedraltypes ] of funct 4
	const std::vector<DihedralBond::Parameters>& GetPeriodicImproperParameters(const std::array<std::string, DihedralBond::nAtoms>& query);

	// The CMAP grid of 5 atomtypes, as an index into GetCmapGrids. Matched in the given order only, as GROMACS does
	std::optional<int> GetCmapGridIndex(const std::array<std::string, CmapBond::nAtoms>& atomtypes) const;
	const std::vector<CmapGrid>& GetCmapGrids() const { return cmapGrids; }

private:
	std::unordered_map<std::string, int> cmapGridIndices;	// Keyed by the 5 atomtypes joined by spaces
	std::vector<CmapGrid> cmapGrids;

	std::unique_ptr<AtomtypeDatabase> ljParameters;

	std::unique_ptr<ParameterDatabase<SinglebondType>> singlebondParameters;
	std::unique_ptr<ParameterDatabase<PairbondType>> pairbondParameters;
	std::unique_ptr<ParameterDatabase<AnglebondType>> anglebondParameters;
	std::unique_ptr<ParameterDatabase<DihedralbondType>> dihedralbondParameters;
	std::unique_ptr<ParameterDatabase<DihedralbondType>> periodicImproperParameters;
	std::unique_ptr<ParameterDatabase<ImproperDihedralbondType>> improperdihedralbondParameters;

	template<typename GenericBond>
	const std::vector<typename GenericBond::Parameters>& _GetBondParameters(const auto& query);

	void LoadFileIntoForcefield(const GenericItpFile& file);
};