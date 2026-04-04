#pragma once

#include "TestUtils.h"
#include <random>

#include "Engine.cuh"

using namespace TestUtils;

namespace KernelAlgorithms {

	LimaUnittestResult WarpSort64_Unittest(EnvMode envmode) {
		bool success = Engine::TestAlgorithms();

		if (!success) {
			return LimaUnittestResult{ false, "Warp sort failed for some bin configs", true };
		}
		return LimaUnittestResult{ true, "Warp sort correct for all bin configs", true };
	}


}