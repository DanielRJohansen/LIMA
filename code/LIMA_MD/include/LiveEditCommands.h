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
	};;

	using Command = std::variant<Invalid, InsertMolecule, DragMolecule>;


	Command ParseCommand(const std::string& commandStr);
}