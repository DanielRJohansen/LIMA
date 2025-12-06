#include <iostream>
#include <filesystem>

#include "CommandlineUtils.h"
#include "MoleculeUtils.h"
#include "Display.h"
#include "argparser.h"


namespace fs = std::filesystem;



int render(int argc, char** argv) {
    const std::string helpText = R"(
Usage: render [OPTIONS]

Description:
    This program renders a .gro file

Options:
    -conf [path]
        Path to the configuration (.gro) file. Defaults to ./conf.gro.
    
    -topology [path]
        Path to the topology (.top) file. Defaults to ./topol.top.

    -whole
        Shows the molecule as awhole and centered, even if the molecule is fragmented due to periodic boundary condition.
        This requires the toplogy file to be provided. Defaults to false.
    
    -hidewater
        Hides water molecules from the rendering. Defaults to false.

    -highlight, -hl [list of idxs]
        Highlights the atoms with the given indices in the rendering. Specify indices, NOT .gro ids.
        Indices are 0-indexed, so the first atom is 0, the second is 1, etc.

Example:
    render -conf myconf.gro -topology mytopol.top -whole -highlight 0 1 4 5
    )";


    ArgParser parser(helpText);

	fs::path conf = "./conf.gro";
	fs::path topol = "./topol.top";
    bool whole = false;
    bool hidewater = false;
    std::vector<int> highlightAtomsInput{};

    parser.AddOption({ "-conf", "-c" }, false, conf);
    parser.AddOption({ "-topology", "-top", "-t" }, false, topol);    
	parser.AddFlag({ "-whole", "-w" }, [&whole]() { whole = true; });
	parser.AddFlag({ "-hidewater", "-hw" }, [&hidewater]() { hidewater = true; });
	parser.AddOption({ "-highlight", "-hl" }, false, highlightAtomsInput);
    parser.Parse(argc, argv);

    GroFile grofile{ conf };

    if (hidewater) {
        while (grofile.atoms.back().residueName == "SOL") {
            grofile.atoms.pop_back();
        }
    }

    if (whole) {
        TopologyFile topfile{ topol };
        MoleculeUtils::CenterMolecule(grofile, topfile.GetMoleculeType());
    }

    std::set<int> highlightedAtoms(highlightAtomsInput.begin(), highlightAtomsInput.end());

    Display d{};
    auto renderTask = std::make_unique<Rendering::GrofileTask>(grofile);
	renderTask->highlightedAtoms = highlightedAtoms;
    d.Render(std::move(renderTask), true);

    return 0;
}
