#include <iostream>
#include <filesystem>
#include <algorithm>
#include <cctype>
#include <format>
#include <memory>
#include <random>
#include <string>
#include <system_error>

#include "CommandlineUtils.h"
#include "MoleculeUtils.h"
#include "Display.h"
#include "argparser.h"
#include "Programs.h"


namespace fs = std::filesystem;

namespace RenderCli {
	class TemporaryDirectory {
	public:
		TemporaryDirectory() {
			const fs::path root = fs::temp_directory_path();
			std::random_device random;
			for (int attempt = 0; attempt < 32; ++attempt) {
				path_ = root / ("lima-render-" + std::to_string(random()) + "-" + std::to_string(attempt));
				std::error_code error;
				if (fs::create_directory(path_, error)) return;
			}
			throw std::runtime_error("Could not create a temporary directory for structure conversion");
		}

		~TemporaryDirectory() {
			std::error_code error;
			fs::remove_all(path_, error);
		}

		TemporaryDirectory(const TemporaryDirectory&) = delete;
		TemporaryDirectory& operator=(const TemporaryDirectory&) = delete;
		const fs::path& path() const { return path_; }

	private:
		fs::path path_;
	};

	inline std::string Lowercase(std::string value) {
		std::ranges::transform(value, value.begin(), [](const unsigned char c) {
			return static_cast<char>(std::tolower(c));
		});
		return value;
	}
}


int render(int argc, char** argv) {
    const std::string helpText = R"(
Usage: lima render [OPTIONS]

Description:
    Render a .gro file, or convert and render a protein .pdb/.cif file.

Options:
    -conf, -f, -structure [path]
        Path to a .gro, .pdb, or .cif file. Defaults to ./conf.gro.
    
    -topology [path]
        Topology for a .gro input. Defaults to ./topol.top. PDB/CIF inputs use
        the topology produced by their internal CHARMM27 conversion.

    -water, -water-model [tip3p|tip4p|tips3p|tip5p|spc|spce]
        Water topology used during PDB/CIF conversion. Default: tip3p.

    -whole
        Shows the molecule as awhole and centered, even if the molecule is fragmented due to periodic boundary condition.
        This requires the toplogy file to be provided. Defaults to false.
    
    -hidewater
        Hides water molecules from the rendering. Defaults to false.

    -highlight, -hl [list of idxs]
        Highlights the atoms with the given indices in the rendering. Specify indices, NOT .gro ids.
        Indices are 0-indexed, so the first atom is 0, the second is 1, etc.

Example:
    lima render -f myconf.gro -topology mytopol.top -whole
    lima render -f protein.cif -whole -highlight 0 1 4 5
    )";


    ArgParser parser(helpText);

	fs::path conf = "./conf.gro";
	fs::path topol = "./topol.top";
	std::string water = "tip3p";
    bool whole = false;
    bool hidewater = false;
    std::vector<int> highlightAtomsInput{};

    parser.AddOption({ "-conf", "-c", "-f", "-structure" }, false, conf);
    parser.AddOption({ "-topology", "-top", "-t" }, false, topol);    
	parser.AddOption({ "-water", "-water-model" }, false, water);
	parser.AddFlag({ "-whole", "-w" }, [&whole]() { whole = true; });
	parser.AddFlag({ "-hidewater", "-hw" }, [&hidewater]() { hidewater = true; });
	parser.AddOption({ "-highlight", "-hl" }, false, highlightAtomsInput);
    parser.Parse(argc, argv);

	std::unique_ptr<RenderCli::TemporaryDirectory> conversionDirectory;
	bool convertedStructure = false;
	const std::string extension = RenderCli::Lowercase(conf.extension().string());
	if (extension == ".pdb" || extension == ".cif") {
		convertedStructure = true;
		conversionDirectory = std::make_unique<RenderCli::TemporaryDirectory>();
		const Programs::WaterModel waterModel = Programs::ParseWaterModel(water);
		const std::optional<std::string> conversionName{ "render_input" };
		const Programs::GmxConversionResult converted = extension == ".pdb"
			? Programs::pdb2gmx(conf, conversionName, waterModel, conversionDirectory->path())
			: Programs::cif2gmx(conf, conversionName, waterModel, conversionDirectory->path());
		conf = converted.gro;
		topol = converted.topology;
	}
	else if (extension != ".gro") {
		throw std::runtime_error(std::format(
			"lima render expects a .gro, .pdb, or .cif input file, got {}", conf.string()));
	}

    GroFile grofile{ conf };

    if (whole) {
        TopologyFile topfile{ topol };
		MoleculeUtils::MakeMoleculeWholeAfterPBCFragmentation(grofile, topfile);
    }

	if (hidewater) {
		while (!grofile.atoms.empty() && grofile.atoms.back().residueName == "SOL") {
			grofile.atoms.pop_back();
		}
	}

	// Converted coordinate files carry crystallographic cell dimensions, which do
	// not necessarily bound the displayed biological structure. A whole structure
	// can likewise extend beyond its former periodic cell after unwrapping.
	if (convertedStructure || whole) MoleculeUtils::FitMoleculeInBox(grofile);

    std::set<int> highlightedAtoms(highlightAtomsInput.begin(), highlightAtomsInput.end());

    Display d{};
    auto renderTask = std::make_unique<Rendering::GrofileTask>(grofile);
	renderTask->highlightedAtoms = highlightedAtoms;
    d.Render(std::move(renderTask), true);

    return 0;
}
