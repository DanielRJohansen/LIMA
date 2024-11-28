#pragma once

#include "MDFiles.h"

namespace Lipids {
	struct Select {
		Select(const std::string& lipidname, const fs::path& workDir, double percentage);

		const bool userSupplied = false;
		const std::string lipidname;
		const double percentage;
		std::shared_ptr<GroFile> grofile;
		std::shared_ptr<TopologyFile> topfile;
	};
	// Since this is a vector of structs with unique_ptrs, it can never be copied, or resized
	using Selection = std::vector<Select>;

	void OrientLipidhead(GroFile&, Float3 desiredOrientation = Float3{0,0,1});

	void OrganizeLipidIntoCompoundsizedSections(GroFile&, TopologyFile::Moleculetype&);

	void _MakeLipids(bool writeToFile, bool displayEachLipidAndHalt);
};

