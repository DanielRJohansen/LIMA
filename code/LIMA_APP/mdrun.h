#pragma once

#include <iostream>

#ifdef __linux__
#include "autorecompile.h"
#endif

int mdrun(int argc, char** argv)
{
#ifdef __linux__
	// Recompile system to fit simparams
	const int compile_failed = SelfRecompile::autoRecompile();
    if (compile_failed) return 1;
#endif


	// Call the program /opt/LIMA/Applications/mdrun with the same arguments
    const std::string home = getenv("HOME");
    std::string command = home + "/LIMA/applications/mdrun";
    for (int i = 2; i < argc; ++i) {    // argv[1] is mdrun
        command += " ";
        command += argv[i];
    }
    std::cout << "executing cmd: " << command << "\n";
    return system(command.c_str());
}
