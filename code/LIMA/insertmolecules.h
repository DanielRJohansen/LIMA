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
    
    -display, -d
        Render the box as Staticbody Energyminimization is carried out

    -help, -h
        Display this help text and exit.

Example:
    insertmolecules -conf_source mysource.gro -top_source mytopology.top -num_insertions 150 -rr -d
    )";


    ArgParser parser(helpText);

    fs::path confSrcPath{};
    fs::path topSrcPath{};
    fs::path confTgtPath{"./conf.gro"};
    fs::path topTgtPath{"./topol.top"};

    bool rotateRandomly = false;
    int nInsertions{};
    bool display = false;

    parser.AddOption({ "-conf_source", "-cs" }, true, confSrcPath);
    parser.AddOption({ "-top_source", "-ts" }, true, topSrcPath);
    parser.AddOption({ "-conf_target", "-ct" }, false, confTgtPath);
    parser.AddOption({ "-top_target", "-tt" }, false, topTgtPath);
    parser.AddOption({ "-num_insertions", "-n" }, true, nInsertions);
    parser.AddFlag({ "-rotate_randomly", "-rr" }, [&rotateRandomly]() {rotateRandomly = true; });
    parser.AddFlag({ "-display", "-d" }, [&display]() {display = true; });
    parser.Parse(argc, argv);

    if (!fs::exists(confSrcPath)) {printf("Invalid conf src");}//TODO make standard,prettier
    if (!fs::exists(topSrcPath)) {printf("Invalid top src");}
    if (!fs::exists(confTgtPath)) {printf("Invalid conf tgt");}

    GroFile groSrc{ confSrcPath };
    auto topSrc = std::make_shared<TopologyFile>(topSrcPath);
    GroFile groTgt{ confTgtPath };
    TopologyFile topTgt{};
    topTgt.SetSystem(topSrc->GetSystem().title + " " + std::to_string(nInsertions));


    SimulationBuilder::InsertSubmoleculesInSimulation(groTgt, topTgt, groSrc, topSrc, nInsertions, rotateRandomly);
    Programs::StaticbodyEnergyMinimize(groTgt, topTgt, display);
    groTgt.printToFile(confTgtPath);
    topTgt.printToFile(topTgtPath);

    return 0;
}
