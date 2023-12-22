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


    // This is not pretty, but we need to make a system call to the actual program. This is because the mdrun option
    // will have to first recombile the program
    const std::string program = argv[1];
    if (program == "mdrun") mdrun(argc, argv);
    else {
        // Any other lima program can be caught directly by LIMA_SERVER/lima
        std::string command = "~/LIMA/applications/mdrun";
        for (int i = 1; i < argc; ++i) {
            command += " ";
            command += argv[i];
        }

        return system(command.c_str());
    }
    //else std::cout << "lima program: \"" << program << "\" is not recognized\n";


    return 0;
}
