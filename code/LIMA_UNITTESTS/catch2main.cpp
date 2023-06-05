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
		REQUIRE(doPoolCompSolBenchmark(envmode));
	}

	TEST_CASE("TestForceCorrectness::SingleBond") {
		REQUIRE(doSinglebondBenchmark(envmode));
	}

	TEST_CASE("TestForceCorrectness::AngleBond") {
		REQUIRE(doAnglebondBenchmark(envmode));
	}

	TEST_CASE("TestForceCorrectness::DihedralBond") {
		REQUIRE(doDihedralbondBenchmark(envmode));
	}

	TEST_CASE("TestForceCorrectness::Methionine") {
		REQUIRE(doMethionineBenchmark(envmode));
	}

	TEST_CASE("TestForceCorrectness::TenSolvents") {
		REQUIRE(TestUtils::loadAndRunBasicSimulation("TenSolvents", envmode, 0.05, {}, false));
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
		InputSimParams ip{ 10.f, 10 };
		TestUtils::loadAndRunBasicSimulation("SolventBenchmark", envmode, 0.002, {ip}, true);

		//InputSimParams ip{ 100.f, 10 };
		////auto env = TestUtils::basicSetup("TenSolvents", {ip}, envmode);
		//auto env = TestUtils::basicSetup("SolventBenchmark", { ip }, envmode);
		//env->prepareForRun();
		//env->run(true);
		//Environment env{ "aa", envmode };
	}
}

namespace StressTesting {
	TEST_CASE("StressTesting::RepeatPool100x") {
		auto func = []() {doPoolBenchmark(envmode); };
		TestUtils::stressTest(func, 3);
	}
}