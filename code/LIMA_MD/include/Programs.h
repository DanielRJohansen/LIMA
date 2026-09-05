#pragma once

#include "SimulationBuilder.h"

#include <string_view>

class Environment;
class MoleculeHullCollection;
class Simulation;

namespace Programs {
	enum class WaterModel { Tip3p, Tip4p, Tips3p, Tip5p, Spc, Spce };

	struct GmxConversionResult {
		fs::path gro;
		fs::path topology;
		fs::path positionRestraints;
	};

	WaterModel ParseWaterModel(std::string_view name);

	void GetForcefieldParams(const GroFile&, const TopologyFile&, const fs::path& workdir);

	MoleculeHullCollection MakeLipidVesicle(GroFile&, TopologyFile&, Lipids::Selection, float vesicleRadius, 
		Float3 vesicleCenter, std::optional<int> numLipids=std::nullopt);

	void MoveMoleculesUntillNoOverlap(MoleculeHullCollection& mhCol, Float3 boxSize, bool renderProgress);

	/// <summary></summary>
	/// <param name="writePositionsToGrofile">If false, the grofile will not be modified</param>
	/// <param name="mayOverlapEdges">If the box contents may spill over the edge, set this to true.
	/// Then we will first run a pre-EM with boxEdgePotential enabled</param>
	/// <returns>Can be discarded if not needed. Only makes sense to discard if overwriting grofile</returns>
	std::unique_ptr<Simulation> EnergyMinimize(GroFile&, const TopologyFile&,
		bool writePositionsToGrofile, const fs::path& workDir, EnvMode, bool mayOverlapEdges, float emtol=100.f);

	void StaticbodyEnergyMinimize(GroFile&, const TopologyFile&, bool render);

	/// Build a CHARMM27 topology and coordinates from a PDB structure. Output is
	/// written next to pdbfile as conf.gro/topol.top/posre.itp, or with the
	/// supplied basename as <name>.gro/<name>.top/<name>_posre.itp. The selected
	/// water model determines the solvent topology included for later solvation.
	GmxConversionResult pdb2gmx(const fs::path& pdbfile, std::optional<std::string> name = std::nullopt,
		WaterModel waterModel = WaterModel::Tip3p, std::optional<fs::path> outputDirectory = std::nullopt);

	/// Build CHARMM27 topology and coordinates from an mmCIF structure.
	GmxConversionResult cif2gmx(const fs::path& ciffile, std::optional<std::string> name = std::nullopt,
		WaterModel waterModel = WaterModel::Tip3p, std::optional<fs::path> outputDirectory = std::nullopt);
}
