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
    std::string command = "~/LIMA/applications/mdrun";
    for (int i = 1; i < argc; ++i) {
        command += " ";
        command += argv[i];
    }

    return system(command.c_str());
}
