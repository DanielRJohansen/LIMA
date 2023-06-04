// Due to the following define, i belive all includes to test-headers must come after
#define CATCH_CONFIG_MAIN

#include "LIMA_TESTS/src/ForceCorrectness.h"
#include "LIMA_TESTS/src/MDStability.cuh"

#include <Catch2/catch.hpp>	// Don't include windows.h after this line

constexpr auto envmode = Environment::Mode::Headless;

namespace TestForceCorrectness {
	TEST_CASE("TestForceCorrectness::Pool") {
		doPoolBenchmark(envmode);			// Two 1-particle molecules colliding
		SUCCEED();
	}
}

namespace TestMDStability {

	TEST_CASE("TestMDStability::EMRunSolvent", "hey") {
		REQUIRE(loadAndEMAndRunBasicSimulation("SolventBenchmark", envmode, 0.0002));
		SUCCEED();
	}

	TEST_CASE("TestMDStability::ahhh") {
		InputSimParams ip{ 100.f, 1 };
		auto env = TestUtils::basicSetup("SolventBenchmark", {ip}, envmode);
		env.prepareForRun();
		env.run(true);
		FAIL();
		//Environment env{ "aa", envmode };
		SUCCEED();
	}
}
