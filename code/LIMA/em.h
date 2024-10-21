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

    bool render = false;

    parser.AddOption({ "-conf", "-c" }, false, confPath);
    parser.AddOption({ "-top", "-t" }, false, topPath);
    parser.AddFlag({ "-display", "-d" }, [&render](){render=true;});
    parser.Parse(argc, argv);

    GroFile grofile{ confPath };
    TopologyFile topfile{ topPath };

    Programs::EnergyMinimize(grofile, topfile, true, fs::current_path(), render ? Full : ConsoleOnly, false);

    grofile.printToFile();
    topfile.printToFile();

    return 0;
}
