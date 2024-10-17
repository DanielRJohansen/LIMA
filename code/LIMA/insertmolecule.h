#include "argparser.h"
#include <filesystem>
#include "MDFiles.h"
#include "SimulationBuilder.h"


int insertmolecule(int argc, char** argv) {
    namespace fs = std::filesystem;

    const std::string helpText = R"(
Usage: insertmolecule [OPTIONS]

Description:
    This programs inserts a molecule in a target box

Options:
    -conf_source, -cs [path]
        Path to the source .gro file. <Required>

    -top_source, -ts [path]
        Path to the source .top file. <Required>
    
    -conf_target, -ct [path]
        Path to the target .gro file. Defaults to ./conf.gro

    -top_target, -tt [path]
		Path to the target .top file. Defaults to ./topol.top
    
    -position, -p [x y z]
		Target position to insert the molecule. Defaults to the center of conf_target.

    -help, -h
        Display this help text and exit.

Example:
    insertmolecule -conf_source mysource.gro -top_source mytopology.top -position 2.14 9 5.1413
    )";


    ArgParser parser(helpText);

    fs::path confSrcPath{};
    fs::path topSrcPath{};
    fs::path confTgtPath{"./conf.gro"};
    fs::path topTgtPath{"./topol.top"};

    Float3 position{FLT_MAX, FLT_MAX , FLT_MAX };

    parser.AddOption({ "-conf_source", "-cs" }, true, confSrcPath);
    parser.AddOption({ "-top_source", "-ts" }, true, topSrcPath);
    parser.AddOption({ "-conf_target", "-ct" }, false, confTgtPath);
    parser.AddOption({ "-top_target", "-tt" }, false, topTgtPath);
    parser.AddOption({ "-position", "-p" }, false, position);


    GroFile groSrc{ confSrcPath };
    auto topSrc = std::make_shared<TopologyFile>(topSrcPath);
    GroFile groTgt{ confTgtPath };
    TopologyFile topTgt{ topTgtPath };


    if (position.x == FLT_MAX) {
        position = groTgt.box_size / 2.f;
	}

    SimulationBuilder::InsertSubmoleculeInSimulation(groTgt, topTgt, groSrc, topSrc, position);

    return 0;
}