#pragma once

#include <iostream>

#ifdef __linux__
#include "autorecompile.h"
#endif

int buildMembrane(int argc, char** argv)
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

    //const fs::path work_dir = simulations_dir + "/BuildMembraneSmall";
    //Environment env{ work_dir.string(), envmode, false };

    //env.CreateSimulation(7.f);
    //env.createMembrane(do_em);

    return system(command.c_str());
}
