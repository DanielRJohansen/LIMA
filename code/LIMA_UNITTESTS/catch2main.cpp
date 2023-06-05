// Due to the following define, i belive all includes to test-headers must come after
#define CATCH_CONFIG_MAIN

#include "LIMA_TESTS/src/ForceCorrectness.h"
#include "LIMA_TESTS/src/MDStability.cuh"

#include <Catch2/catch.hpp>	// Don't include windows.h after this line

constexpr auto envmode = Environment::Mode::Headless;

namespace TestForceCorrectness {
	TEST_CASE("TestForceCorrectness::Pool") {
		REQUIRE(doPoolBenchmark(envmode));			// Two 1-particle molecules colliding
	}

	TEST_CASE("TestForceCorrectness::PoolCarbonSol") {
		doPoolCompSolBenchmark(envmode);
	}

	TEST_CASE("TestForceCorrectness::SingleBond") {
		doSinglebondBenchmark(envmode);
	}

	TEST_CASE("TestForceCorrectness::AngleBond") {
		doAnglebondBenchmark(envmode);
	}

	TEST_CASE("TestForceCorrectness::DihedralBond") {
		doDihedralbondBenchmark(envmode);
	}

	TEST_CASE("TestForceCorrectness::Methionine") {
		doMethionineBenchmark(envmode);
	}
}

namespace TestMDStability {

	TEST_CASE("TestMDStability::EMRunSolvent", "hey") {
		REQUIRE(loadAndEMAndRunBasicSimulation("SolventBenchmark", envmode, 0.0002));
		SUCCEED();
	}

	TEST_CASE("TestMDStability::EightResiduesNoSolvent", "hey") {
		REQUIRE(doEightResiduesNoSolvent(envmode, 0.0002));		
	}

	TEST_CASE("TestMDStability::ahhh") {
		InputSimParams ip{ 100.f, 10 };
		auto env = TestUtils::basicSetup("TenSolvents", {ip}, envmode);
		//auto env = TestUtils::basicSetup("TorsionBenchmark", { ip }, envmode);
		env->prepareForRun();
		env->run(true);
		//Environment env{ "aa", envmode };
	}
}
