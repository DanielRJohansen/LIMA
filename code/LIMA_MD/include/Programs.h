#pragma once

#include "SimulationBuilder.h"
#include "Environment.h"

#include <string_view>

class MoleculeHullCollection;

namespace Programs {
	enum class WaterModel { Tip3p, Tip4p, Tips3p, Tip5p, Spc, Spce };

	struct GmxConversionResult {
		GroFile grofile;
		TopologyFile topology;
		// Ordered like topology.GetSystem().molecules. Each file contains the
		// position restraints for the corresponding molecule.
		std::vector<GenericItpFile> positionRestraints;
		// What the conversion decided or dropped (disulfides, histidine states, skipped records), for the user
		std::vector<std::string> notes;
	};

	struct ToGmxOptions {
		// AMBER99SB only has charged termini: its neutral and capped ends are separate residues (ACE, NME, NHE)
		enum class ForceField { Charmm27, Amber99sb } forceField = ForceField::Charmm27;
		WaterModel waterModel = WaterModel::Tip3p;
		// Histidine protonation: Hsd/Hse/Hsp for all histidines, or Auto from the H-bond network.
		// Explicit input names (HISD, HISE, HISH) always win
		enum class His { Hsd, Hse, Hsp, Auto } his = His::Hsd;
		enum class NTerminus { Charged, Neutral } nTerminus = NTerminus::Charged;		// NH3+, NH2
		enum class CTerminus { Charged, Neutral, Amide, Methylamide } cTerminus = CTerminus::Charged;	// COO-, COOH, CT2, CT3
		bool disulfides = true;		// Bond CYS pairs whose SG atoms are within bonding distance
		bool ignoreHydrogens = false;	// Drop the input's hydrogens and build them all
		char altloc = 'A';			// Alternate location to keep
		std::vector<std::string> chains;	// Chains to keep, empty for all
	};

	WaterModel ParseWaterModel(std::string_view name);
	// Filename of the model's topology in a forcefield dir, eg. "tip3p.itp"
	std::string_view WaterModelItpFilename(WaterModel model);

	void GetForcefieldParams(const GroFile&, const TopologyFile&, const fs::path& workdir);

	MoleculeHullCollection MakeLipidVesicle(GroFile&, TopologyFile&, Lipids::Selection, float vesicleRadius, 
		Float3 vesicleCenter, std::optional<int> numLipids=std::nullopt);

	void MoveMoleculesUntillNoOverlap(MoleculeHullCollection& mhCol, Float3 boxSize, bool renderProgress);

	void StaticbodyEnergyMinimize(GroFile&, const TopologyFile&, bool render);

	// The lipids are loaded and the membrane built in the job's preprocess, on Environment's preparation workers
	SimulationJob MakeMembraneJob(fs::path workDir, Lipids::Composition composition,
		Float3 boxSize, MembraneGeometry::Figure geometry, int seed,
		SimParams params, EnvMode mode = EnvMode::Headless, bool solvate = false);

	SimulationJob MakeSimulationJob(fs::path workDir, MolecularSystem system,
		SimParams params, EnvMode mode = EnvMode::Headless);

	/// Build in-memory CHARMM27 or AMBER99SB coordinates, topology, and position restraints
	/// from a protein PDB or mmCIF structure.
	GmxConversionResult ToGmx(const fs::path& structureFile, WaterModel waterModel = WaterModel::Tip3p);
	GmxConversionResult ToGmx(const fs::path& structureFile, const ToGmxOptions& options);
}
