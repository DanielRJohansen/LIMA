#pragma once

#include "TestUtils.h"
#include "Environment.h"
#include "Printer.h"
#include "Utilities.h"
#include "LimaTypes.cuh"
//#include "EngineUtils.cuh"



#include <iostream>
#include <string>
#include <algorithm>


namespace ElectrostaticsTests {
	using namespace TestUtils;

	static LimaUnittestResult TestChargedParticlesEndInCorrectSection(EnvMode envmode, float max_vc = 0.05, float max_gradient = 1e-5) {
		SimParams simparams{ 50000, 20, true, PBC };
		simparams.coloring_method = ColoringMethod::Charge;
		simparams.data_logging_interval = 20;
		simparams.snf_select = HorizontalChargeField;
		auto env = basicSetup("ElectrostaticField", { simparams }, envmode);


		// Do sim
		env->run();
		
		return LimaUnittestResult{ LimaUnittestResult::SUCCESS, "", envmode == Full };
	}
}

