#include "argparser.h"

int makebox(int argc, char** argv) {
    const std::string helpText = R"(
Usage: makebox [OPTIONS]

Description:
    This program makes an empty box (.gro & .top file). Can only create cubic boxes, as this is a programwide limitation enforced by lima (for performance)

Options:
    -boxsize [int]
        Size of the box in nm. Defaults to 10 nm. <Required>

    -name [string]
        The name for the output conf and top file. Defaults to conf

    -top_name [path]
        Path to output top file. Defaults to ./topol.top

    -help, -h
        Display this help text and exit.

Example:
    makebox -name mybox -boxsize 14
    )";


    ArgParser parser(helpText);

    std::string name ="";
    int boxsize{};

    // Add options to parser
    parser.AddOption({ "-boxsize", "-b" }, true, boxsize);
    parser.AddOption({ "-name", "-n"}, false, name);
    
    parser.Parse(argc, argv);

	const fs::path confPath = fs::absolute(!name.empty() ? fs::path(name + ".gro") : "conf.gro");
	const fs::path topPath  = fs::absolute(!name.empty() ? fs::path(name + ".top") : "topol.top");

    //printf("Creating files: \n\t%s\n\t%s\n", confPath.string(), topPath.string());


    // Use the paths in your program
    GroFile grofile{};
    grofile.m_path = confPath;
    grofile.box_size = Float3{ static_cast<float>(boxsize) };
    grofile.printToFile();

    TopologyFile topfile{};
    topfile.path = topPath;
    topfile.printToFile();

    return 0;
}