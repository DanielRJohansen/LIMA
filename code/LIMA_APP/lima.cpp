#include <iostream>

#include "mdrun.h"


void sayHello() {
    std::printf("Welcome to LIMA Dynamics\n");
}

int main(int argc, char** argv) {

    sayHello();


    if (argc == 1) {
        std::printf("Please call lima with a program to execute (e.g. lima mdrun conf.gro topol.top)");
        return 0;
    }

    const std::string program = argv[1];
    if (program == "mdrun") mdrun(argc, argv);
    else std::cout << "lima program: \"" << program << "\" is not recognized\n";


    return 0;
}
