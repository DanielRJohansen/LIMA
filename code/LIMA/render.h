#include <iostream>
#include <filesystem>

#include "CommandlineUtils.h"
#include "MoleculeUtils.h"
#include "Display.h"

namespace fs = std::filesystem;

struct RenderSetup {
    RenderSetup(int argc, char** argv) {
        for (int i = 1; i < argc; ++i) {
            std::string arg = CmdLineUtils::ToLowercase(argv[i]);

            if (arg == "-conf") {
                if (i + 1 < argc) {
                    conf = argv[++i];
                }
                else {
                    std::cerr << "-conf expects a path argument." << std::endl;
                    exit(1);
                }
            }
            else if (arg == "-topology") {
                if (i + 1 < argc) {
                    topol = argv[++i];
                }
                else {
                    std::cerr << "-topology expects a path argument." << std::endl;
                    exit(1);
                }
            }
            else if (arg == "-whole")
                whole = true;                
            else if (arg == "-hidewater")
                hideWater = true;
            else if (arg == "-help" || arg == "-h") {
                std::cout << helpText;
                exit(0);
            }
            else {
                std::cerr << "Unknown argument: " << arg << std::endl;
                exit(1);
            }
        }
    }

    fs::path conf = "./conf.gro";
    fs::path topol = "./topol.top";
    bool whole = false;
    bool hideWater = true;
private:
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

Example:
    render -conf myconf.gro -topology mytopol.top -whole
    )";
};

int render(int argc, char** argv) {

    RenderSetup setup{ argc, argv };

    GroFile grofile{ setup.conf };

    if (setup.hideWater) {
        while (grofile.atoms.back().residueName == "SOL") {
            grofile.atoms.pop_back();
        }
    }

    if (setup.whole) {
        TopologyFile topfile{ setup.topol };
        MoleculeUtils::CenterMolecule(grofile, topfile.GetMoleculeType());
    }
    Display d{ Full };
    d.Render(std::make_unique<Rendering::GrofileTask>(grofile), true);

    return 0;
}
