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
        Path to output gro file. Defaults to ./solvated.gro
        
    -pressure [int]
        Defaults to 34 solvents/nm3. 

    -help, -h
        Display this help text and exit.

Example:
    lima solvate -conf conf.gro
    )";


    ArgParser parser(helpText);

    fs::path input = "./conf.gro";
	fs::path output = "./solvated.gro";
	int pressure = SimulationBuilder::defaultSolventsPerNm3; // [solvents/nm3]

    // Add options to parser
    parser.AddOption({ "-conf", "-c", "-input"}, false, input);
	parser.AddOption({ "-output", "-o", "-outputconf" }, false, output);
    parser.AddOption({ "-pressure", "-p"}, false, pressure);
    parser.Parse(argc, argv);

    if (!fs::exists(input)) {
        std::cerr << "Input file does not exist: " << input.string() << "\n";
        return 1;
	}

    // Use the paths in your program
    GroFile grofile{ input };

    SimulationBuilder::SolvateGrofile(grofile);

	grofile.printToFile(output);

    return 0;
}