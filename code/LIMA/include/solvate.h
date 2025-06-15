#include "argparser.h"
#include "SimulationBuilder.h"

int solvate(int argc, char** argv) {
    const std::string helpText = R"(
Usage: makebox [OPTIONS]

Description:
    This program makes an empty box (.gro & .top file). Can only create cubic boxes, as this is a programwide limitation enforced by lima (for performance)

Options:

    -conf [.gro file]
        Path to input conf file. Defaults to ./conf.gro
    
    -output [.gro file]
        Path to output gro file. Defaults to ./conf_solvated.gro
        
    -pressure [??/??]
        Not yet implemented

    -help, -h
        Display this help text and exit.

Example:
    lima solvate -conf conf.gro
    )";


    ArgParser parser(helpText);

    fs::path input = "./conf.gro";
	fs::path output = "./conf_solvated.gro";

	float pressure = 0.f;  // Placeholder, not yet implemented

    // Add options to parser
    parser.AddOption({ "-conf", "-c", "input"}, false, input);
	parser.AddOption({ "-output", "-o", "outputconf" }, false, output);
    parser.AddOption({ "-pressure", "-p"}, false, pressure);
    parser.Parse(argc, argv);

    // Use the paths in your program
    GroFile grofile{input};

    SimulationBuilder::SolvateGrofile(grofile); // TODO: This only generates O, needs the H2...

	grofile.printToFile(output);

    return 0;
}