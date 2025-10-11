#include "argparser.h"
#include "SimulationBuilder.h"

int solvate(int argc, char** argv) {
    const std::string helpText = R"(
Usage: makebox [OPTIONS]

Description:
    This program makes an empty box (.gro & .top file). Can only create cubic boxes, as this is a programwide limitation enforced by lima (for performance)

Options:

    -confInput, -ci [.gro file]
        Path to input configuration file. Defaults to ./conf.gro
    
    -topInput, -ti [.top file]
        Path to input topology file. Defaults to ./topol.top
        
    -confOutput, -co [.gro file]
        Path to output configuration file. Defaults to {confInput}_solvated.gro

    -topOutput, -to [.top file]
        Path to output topology file. Defaults to {topInput}_solvated.top
    -pressure [int]
        Defaults to 34 solvents/nm3. 

    -help, -h
        Display this help text and exit.

Example:
    lima solvate -conf conf.gro
    )";


    ArgParser parser(helpText);

    fs::path inputConf = "./conf.gro";
    fs::path outputConf{};
    fs::path inputTop = "topol.top";
    fs::path outputTop{};
	int pressure = SimulationBuilder::defaultSolventsPerNm3; // [solvents/nm3]

    // Add options to parser
    parser.AddOption({ "-confInput", "-ci"}, false, inputConf);
	parser.AddOption({ "-topInput", "-ti"}, false, inputTop);
	parser.AddOption({ "-confOutput", "-co" }, false, outputConf);
	parser.AddOption({ "-topOutput", "-to" }, false, outputTop);
    parser.AddOption({ "-pressure", "-p"}, false, pressure);
    parser.Parse(argc, argv);

    if (!fs::exists(inputConf)) {
        std::cerr << "Input file does not exist: " << inputConf.string() << "\n";
        return 1;
	}    

    if (outputConf.empty())
		outputConf = inputConf.stem().string() + "_solvated.gro";
	if (outputTop.empty())
		outputTop = inputTop.stem().string() + "_solvated.top";

    TopologyFile topfile = fs::exists(inputTop) 
		? TopologyFile{ inputTop }
	    : TopologyFile{};
        
    if (!fs::exists(inputTop))
		topfile.title = "Solvated system";
    topfile.path = outputTop;

    // Use the paths in your program
    GroFile grofile{ inputConf };
    //TopologyFile topfile{};    
    SimulationBuilder::SolvateGrofile(grofile, topfile);

	grofile.printToFile(outputConf);
	topfile.printToFile(outputTop);

    return 0;
}