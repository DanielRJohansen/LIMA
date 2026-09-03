#pragma once

#include "LimaTypes.cuh"
#include "MembraneGeometry.h"
#include <filesystem>
#include <optional>
#include <variant>

namespace LiveEdit {
	namespace fs = std::filesystem;

	struct Invalid {};

	struct InsertMolecule {
		fs::path groPath;
		fs::path topPath;
		std::optional<Float3> position = std::nullopt;
	};

	struct MoveMolecule { // TODO: Rename to MoveSelection
		Float3 draggingForce{};
		Float3 rotation{};
	};

	struct BuildMembrane {
		std::vector<std::tuple<std::string, double>> lipids; // {name, percentage}
		std::optional<MembraneGeometry::Figure> geometry = std::nullopt;
	};

	struct TogglePause {};
	struct StepOnce {};	

	struct AtomSelected {
		int particleId = -1;// Global id 
	};

	struct SelectAtomsBasedOnQualifier {
		enum class Qualifier {
			All, Solvent, Nonsolvent
		};		
		Qualifier qualifier;
	};

	struct AddForcemaskToSelection {
		Float3 forcemask{};
	};

	struct ElasticPosition {
		bool x = false;
		bool y = false;
		bool z = false;
	};

	struct EnergyMinimize {
		//float targetMaxforce = 1e+3; // [kJ/mol/nm]
	};

	using Command = std::variant<Invalid, InsertMolecule, MoveMolecule, BuildMembrane,
		TogglePause, AtomSelected, SelectAtomsBasedOnQualifier, AddForcemaskToSelection, ElasticPosition, EnergyMinimize, StepOnce>;


	Command ParseCommand(const std::string& commandStr);
}
