#include "argparser.h"

int makebox(int argc, char** argv) {
    const std::string helpText = R"(
Usage: makebox [OPTIONS]

Description:
    This program makes an empty box (.gro & .top file). Can only create cubic boxes, as this is a programwide limitation enforced by lima (for performance)

Options:
    -boxsize [int]
        Size of the box in nm. Defaults to 10 nm. <Required>

    -conf_name [path]
        Path to output conf file. Defaults to ./conf.gro

    -top_name [path]
        Path to output top file. Defaults to ./topol.top

    -help, -h
        Display this help text and exit.

Example:
    makebox -name mybox -boxsize 14
    )";


    ArgParser parser(helpText);

    fs::path confname = "./conf.gro";
    fs::path topname = "./topol.top";
    int boxsize{};

    // Add options to parser
    parser.AddOption({ "-boxsize" }, true, boxsize);
    parser.AddOption({ "-conf_name" }, false, confname);
    parser.AddOption({ "-top_name" }, false, topname);
    parser.Parse(argc, argv);

    // Use the paths in your program
    GroFile grofile{};
    grofile.m_path = fs::absolute(confname);  // Convert to absolute path
    grofile.box_size = Float3{ static_cast<float>(boxsize) };
    grofile.printToFile();

    TopologyFile topfile{};
    topfile.path = fs::absolute(topname);  // Convert to absolute path
    topfile.printToFile();

    return 0;
}