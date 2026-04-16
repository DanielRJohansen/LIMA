#pragma once

#include "LimaTypes.cuh"
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

	struct DragMolecule {
		int particleId = -1;// Global id of any particle in the molecule
		Float3 draggingForce{};
	};

	struct BuildMembrane {
		std::vector<std::tuple<std::string, double>> lipids; // {name, percentage}
		std::optional<float> membraneCenterZ = std::nullopt;
	};

	struct TogglePause {};
	
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

	using Command = std::variant<Invalid, InsertMolecule, DragMolecule, BuildMembrane, 
		TogglePause, AtomSelected, SelectAtomsBasedOnQualifier, AddForcemaskToSelection>;


	Command ParseCommand(const std::string& commandStr);
}