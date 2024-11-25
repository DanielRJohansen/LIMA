#pragma once

#include <string.h>
#include <fstream>
#include <vector>
#include <array>
#include <memory>

#include "Bodies.cuh"
#include "Constants.h"
#include "Utilities.h"
#include "BoundaryConditionPublic.h"
#include "EngineCore.h"

#include "MDFiles.h"

#include "tuple"
#include <set>
struct BoxImage;





// --------------------------------- Bond Factories --------------------------------- //

struct ParticleFactory {
	const int activeLjtypeParameterIndex = -1;
	const TopologyFile::AtomsEntry* topAtom = nullptr;
	//const GroRecord* groAtom = nullptr;
	const Float3 position;
	const int uniqueResId;
	const int indexInGrofile; // 0-indexed
};

template <int n_Atoms, typename ParamsType>
struct BondFactory {
	static const int nAtoms = n_Atoms;
	BondFactory(const std::array<uint32_t, nAtoms>& ids, const ParamsType& parameters) 
		: params(parameters), global_atom_indexes(ids) {}

	ParamsType params;
	std::array<uint32_t, nAtoms> global_atom_indexes;
	std::string sourceLine;
};
using SingleBondFactory = BondFactory<2, SingleBond::Parameters>;
using AngleBondFactory = BondFactory<3, AngleUreyBradleyBond::Parameters>;
using DihedralBondFactory = BondFactory<4, DihedralBond::Parameters>;
using ImproperDihedralBondFactory = BondFactory<4, ImproperDihedralBond::Parameters>;

struct ParticleToCompoundMapping {
	int compoundId = -1;
	int localIdInCompound = -1;
};
struct ParticleToBridgeMapping {
	int bridgeId = -1;
	int localIdInBridge = -1;
};
using ParticleToCompoundMap = std::vector<ParticleToCompoundMapping>;
using ParticleToBridgeMap = std::vector<std::optional<ParticleToBridgeMapping>>;

namespace LIMA_MOLECULEBUILD {
	struct Topology {
		std::vector<ParticleFactory> particles;

		std::vector<SingleBondFactory> singlebonds;
		std::vector<AngleBondFactory> anglebonds;
		std::vector<DihedralBondFactory> dihedralbonds;
		std::vector<ImproperDihedralBondFactory> improperdihedralbonds;
	};



	std::unique_ptr<BoxImage> buildMolecules(
		const GroFile& gro_file,
		const TopologyFile& top_file,
		VerbosityLevel vl,
		std::unique_ptr<LimaLogger>,
		bool ignore_hydrogens,
		const SimParams& simparams
	);
}

class CompoundFactory : public Compound, public CompoundInterimState {
public:
	CompoundFactory() {}
	CompoundFactory(const int id) : 
		id(id)
	{
		//memset(forceEnergyFarneighborShortrange, 0, sizeof(forceEnergyFarneighborShortrange));
		//memset(forceEnergyImmediateneighborShortrange, 0, sizeof(forceEnergyImmediateneighborShortrange));
		//memset(forceEnergyBonds, 0, sizeof(forceEnergyBonds));
		//memset(forceEnergyBridge, 0, sizeof(forceEnergyBridge));
		/*memset(forces_interim, 0, sizeof(forces_interim));
		memset(potE_interim, 0, sizeof(potE_interim));*/
	}

	void addParticle(const ParticleFactory&,int global_id, float boxlen_nm, BoundaryConditionSelect bc);

	bool hasRoomForRes(int n_particles_in_res) const {					// TODO: Implement, that it checks n atoms in res
		return ((int)n_particles + n_particles_in_res) <= MAX_COMPOUND_PARTICLES;
	}

	void addIdOfBondedCompound(int id);

	void AddBondgroupReference(int particleId, const BondgroupRef& bgRef);

	int id = -1;	// unique lima id

	Float3 positions[MAX_COMPOUND_PARTICLES];	// Extern positions [nm]
	int global_ids[MAX_COMPOUND_PARTICLES]{};		// For debug ddont like this TODO TODO DELETE

	int indicesInGrofile[MAX_COMPOUND_PARTICLES];	// Temp prolly, used to map compounds atoms back to their index in grofile
};


struct TinyMolFactory {
	TinyMolFactory() {}
    TinyMolFactory(const Float3& pos, int tinymolTypeIndex, const std::string& atomType="", Float3 velocity = Float3{})
        : position(pos), state(TinyMolState{ state.vel_prev,Float3{},tinymolTypeIndex }), atomType(atomType)
	{}
	Float3 position;
	TinyMolState state;
    std::string atomType; // Debug only
};


class BondGroupFactory : public BondGroup {

	
	void AddBondParticles(const ParticleToCompoundMap&, std::span<const uint32_t> globalIds);
public:
	BondGroupFactory() {}

	bool HasSpaceForParticlesInBond(const std::span<const uint32_t>& particleIds) const;

	//void AddParticles(const std::span<const uint32_t>& particleIds);

	void AddBond(const ParticleToCompoundMap&, const SingleBondFactory&);
	void AddBond(const ParticleToCompoundMap&, const AngleBondFactory&);
	void AddBond(const ParticleToCompoundMap&, const DihedralBondFactory&);
	void AddBond(const ParticleToCompoundMap&, const ImproperDihedralBondFactory&);
	
	std::array<uint32_t, maxParticles> particleGlobalIds;
	std::unordered_map<uint32_t, uint8_t> particleGlobalToLocalId;

	static std::vector<BondGroupFactory> MakeBondgroups(const LIMA_MOLECULEBUILD::Topology&,
		const std::vector<ParticleToCompoundMapping>& particlesToCompoundIdMap);

	static std::vector<std::set<BondgroupRef>> MakeParticleToBondgroupsMap(
		const std::vector<BondGroupFactory>&, int nParticlesTotal);

	static std::vector<BondGroup> FinishBondgroups(const std::vector<BondGroupFactory>&);
};



// A translation unit between Gro file representation, and LIMA Box representation
struct BoxImage {
	const std::vector<CompoundFactory> compounds;
	const int32_t total_compound_particles;

	std::vector<BondedParticlesLUT> bpLutCollection;
	
	const std::vector<TinyMolFactory> solvent_positions;

	GroFile grofile;

	const ForceField_NB forcefield;

	const ForcefieldTinymol tinymolTypes;

	LIMA_MOLECULEBUILD::Topology topology; // This is only used for debugging purposes

	const std::vector<NonbondedInteractionParams> nonbondedInteractionParams;

	const std::vector<BondGroup> bondgroups;
};
 
