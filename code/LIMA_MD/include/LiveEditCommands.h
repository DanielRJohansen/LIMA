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

	

	using Command = std::variant<Invalid, InsertMolecule>;


	Command ParseCommand(const std::string& commandStr);
}