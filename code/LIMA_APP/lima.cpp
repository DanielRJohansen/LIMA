// This is basically a wrapper for the real lima program
// Here it is determined if we need to recompile LIMA, and do the recompile if needed. 
// Finally we simply forward all args to the limaserver program
#include <stdlib.h>
#include <iostream>
#include <filesystem>
#include "autorecompile.h"

namespace fs = std::filesystem;
const std::string limaserver_path = std::string(getenv("HOME")) + "/LIMA/applications/limaserver";


constexpr bool requiresRecompile(const std::string& program)
{
    // The first time we call LIMA, we wont have the limaserver yet,
    // so yes recompile
    if (!fs::exists(limaserver_path))
        return true;


    const bool requiresEngine = program == "mdrun" || program == "buildmembrane";

    // TODO: Add some logic to check if the users params are the same as the current engine is compiled with
    return requiresEngine;
}

int main(int argc, char** argv) 
{
    if (argc < 2) {
        std::printf("Please call lima with a program to execute (e.g. lima mdrun conf.gro topol.top)\n");
        return 0;
    }

    const std::string program = argv[1];

    if (requiresRecompile(program)) {
        const int compile_failed = SelfRecompile::autoRecompile();
        if (compile_failed) return 1;
    }

    std::string command = limaserver_path; // Use getenv(home insteadmv
    for (int i = 1; i < argc; ++i) {
        command += " ";
        command += argv[i];
    }

    return system(command.c_str());
}
