#pragma once

#include <iostream>

#ifdef __linux__
#include "autorecompile.h"
#endif

int mdrun(int argc, char** argv) {

#ifdef __linux__
	// Recompile system to fit simparams
	SelfRecompile::autoRecompile();
#endif

	std::printf("mdrun finished\n");
	return 0;
}
