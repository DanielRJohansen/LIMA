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

	std::string NameSelection(const Selection& selection);

	// A lipid mixture by name. Unlike a Selection it reads no files, so it is cheap to create anywhere and is loaded
	// where the membrane is built, eg. in a SimulationJob::preprocess on Environment's preparation workers
	struct Composition {
		struct Lipid { std::string name; double percentage; };
		std::vector<Lipid> lipids;
		fs::path workDir; // <name>.gro/.itp found here take precedence over the stock lipids

		Selection Load() const;
		std::string Name() const; // Same as NameSelection of the loaded Selection
	};

	void OrientLipidhead(GroFile&, Float3 desiredOrientation = Float3{0,0,1});

	void OrganizeLipidIntoCompoundsizedSections(GroFile&, TopologyFile::Moleculetype&);

	// These functions are not for runtime, they prepare data once and for all
	void _MakeLipid(const std::string& name);
	void _MakeLipids(bool writeToFile, bool displayEachLipidAndHalt);

};

