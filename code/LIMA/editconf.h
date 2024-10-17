#include "argparser.h"
#include <filesystem>
#include "MDFiles.h"
#include "SimulationBuilder.h"


int editconf(int argc, char** argv) {
    const std::string helpText = R"(
Usage: insertmolecule [OPTIONS]

Description:
    This programs edits a molecule. Undefined behaviour if the .gro file contains multiple molecules. 
    Transformations are applied in the order in which they apear in this text.


Options:
    -conf, -c [path]
        Path to the input .gro file. <Required>

    -top, -t [path]
        Path to the input .top file. <Required>. This file will not be modified.
    
    -conf_out, -co [path]
        Path to the output .gro file. Defaults to same as input

    -whole, -w
        Will make a molecule that was fragmented due to PBC whole again.
    
    -set_center, -sc [x y z]
		Set the geometric center of the molecule.

    -rotate, -r [x y z]
		Rotate the molecule around the z-axis, y-axis and finally x-axis. Angles in radians.

    -help, -h
        Display this help text and exit.

Example:
    //TODO make example
    )";


    ArgParser parser(helpText);

    fs::path confInputPath{};
    fs::path topInputPath{};
    fs::path confOutputPath{};    

    bool makeWhole = false;
    Float3 setCenter{FLT_MAX, FLT_MAX , FLT_MAX };
    Float3 rotate{FLT_MAX, FLT_MAX , FLT_MAX };

    parser.AddOption({ "-conf", "-c" }, true, confInputPath);
    parser.AddOption({ "-top", "-t" }, true, topInputPath);
    parser.AddOption({ "-conf_out", "-co" }, false, confOutputPath);
    parser.AddFlag({ "-whole", "-w" }, [&makeWhole]() {makeWhole = true; });
    parser.AddOption({ "-set_center", "-sc" }, false, setCenter);
    parser.AddOption({ "-rotate", "-r" }, false, rotate);
    parser.Parse(argc, argv);

    if (confOutputPath.empty())
        confOutputPath = confInputPath;
    

    GroFile grofile{ confInputPath };    
    TopologyFile topfile{ topInputPath };

    bool setCenterChosen = setCenter.x != FLT_MAX;
    bool rotateChosen = rotate.x != FLT_MAX;

    if (makeWhole || setCenterChosen || rotateChosen) {
        MoleculeUtils::MakeMoleculeWholeAfterPBCFragmentation(grofile, topfile);
    }

    if (setCenterChosen) {
        MoleculeUtils::CenterMolecule(grofile, topfile, setCenter);
    }

    if (rotateChosen) {
        MoleculeUtils::RotateMolecule(grofile, rotate);
	}

    grofile.printToFile(confOutputPath);

    return 0;
}