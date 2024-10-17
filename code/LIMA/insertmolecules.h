#include "argparser.h"
#include <filesystem>
#include "MDFiles.h"
#include "SimulationBuilder.h"


int insertmolecules(int argc, char** argv) {
    namespace fs = std::filesystem;

    const std::string helpText = R"(
Usage: insertmolecule [OPTIONS]

Description:
    This programs inserts a molecule multiple times in a target box, with a uniform distribution.
    To avoid gross overlap, the program performs a short static-body energy-minimization.

Options:
    -conf_source, -cs [path]
        Path to the source .gro file. <Required>

    -top_source, -ts [path]
        Path to the source .top file. <Required>
    
    -conf_target, -ct [path]
        Path to the target .gro file. Defaults to ./conf.gro

    -top_target, -tt [path]
		Path to the target .top file. Defaults to ./topol.top

    -num_insertions, -n [int]
        Number of times to insert the molecule. <Required>

    -rotate_randomly, -rr 
		Apply a random rotation to each molecule.

    -help, -h
        Display this help text and exit.

Example:
    insertmolecules -conf_source mysource.gro -top_source mytopology.top 
    )";


    ArgParser parser(helpText);

    fs::path confSrcPath{};
    fs::path topSrcPath{};
    fs::path confTgtPath{"./conf.gro"};
    fs::path topTgtPath{"./topol.top"};

    bool rotateRandomly = false;
    int nInsertions{};

    parser.AddOption({ "-conf_source", "-cs" }, true, confSrcPath);
    parser.AddOption({ "-top_source", "-ts" }, true, topSrcPath);
    parser.AddOption({ "-conf_target", "-ct" }, false, confTgtPath);
    parser.AddOption({ "-top_target", "-tt" }, false, topTgtPath);
    parser.AddOption({ "-num_insertions", "-n" }, true, nInsertions);
    parser.AddFlag({ "-rotate_randomly", "-rr" }, [&rotateRandomly]() {rotateRandomly = true; });


    GroFile groSrc{ confSrcPath };
    auto topSrc = std::make_shared<TopologyFile>(topSrcPath);
    GroFile groTgt{ confTgtPath };
    TopologyFile topTgt{ topTgtPath };

    SimulationBuilder::InsertSubmoleculesInSimulation(groTgt, topTgt, groSrc, topSrc, nInsertions, rotateRandomly);

    // TODO: Implement the staticbody EM

    return 0;
}