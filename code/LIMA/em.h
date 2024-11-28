#include "argparser.h"
#include <filesystem>
#include "MDFiles.h"
#include "SimulationBuilder.h"


int em(int argc, char** argv) {
    namespace fs = std::filesystem;

    const std::string helpText = R"(
Usage: em [OPTIONS]

Description:
    EnergyMinimize a simulation with standard parameters. For more control of the EM, use instead a normal mdrun
    with energyminimize set in the sim_params.txt

Options:
    -conf, -c [path]
        Path to the source .gro file. Defaults to ./conf.gro

    -top, -t [path]
        Path to the source .top file. Defaults to ./topol.top

    -conf_output, -co [path]
        Path to the output .gro file. Defaults to ./conf.gro

    -emtol [value]
        Sets the force tolerance (kJ/mol/nm) for the energy minimization. Default is 100.

    -display -d
        Render the simulation during em

    -help, -h
        Display this help text and exit.

Example:
    em -d
    )";


    ArgParser parser(helpText);

    fs::path confPath{"conf.gro"};
    fs::path topPath{"topol.top"};
    fs::path confPathOut{"conf.gro"};

    bool render = false;
    float emtol = 100.f;

    parser.AddOption({ "-conf", "-c" }, false, confPath);
    parser.AddOption({ "-top", "-t" }, false, topPath);
    parser.AddOption({ "-conf_output", "-co" }, false, confPathOut);
    parser.AddOption({"-emtol"}, false, emtol);
    parser.AddFlag({ "-display", "-d" }, [&render](){render=true;});
    parser.Parse(argc, argv);

    GroFile grofile{ fs::canonical(confPath) };
    TopologyFile topfile{ fs::canonical(topPath) };

    Programs::EnergyMinimize(grofile, topfile, true, fs::current_path(), render ? Full : ConsoleOnly, false, emtol);

    grofile.printToFile(fs::absolute(confPathOut));

    return 0;
}
