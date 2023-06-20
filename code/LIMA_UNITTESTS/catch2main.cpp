// Due to the following define, i belive all includes to test-headers must come after
#define CATCH_CONFIG_MAIN

#include "LIMA_TESTS/src/ForceCorrectness.h"
#include "LIMA_TESTS/src/MDStability.cuh"

#include <Catch2/catch.hpp>	// Don't include windows.h after this line

constexpr auto envmode = EnvMode::Headless;

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

	TEST_CASE("TestForceCorrectness::ImproperDihedralBond") {
		REQUIRE(doImproperDihedralBenchmark(envmode));
	}

	TEST_CASE("TestForceCorrectness::Methionine") {
		REQUIRE(doMethionineBenchmark(envmode));
	}

	TEST_CASE("TestForceCorrectness::TenSolvents") {
		REQUIRE(TestUtils::loadAndRunBasicSimulation("TenSolvents", envmode, 0.0005f));
	}
}

namespace TestMDStability {

	TEST_CASE("TestMDStability::EMRunSolvent") {
		REQUIRE(loadAndEMAndRunBasicSimulation("SolventBenchmark", envmode, 0.0002f));
	}

	TEST_CASE("TestMDStability::EightResiduesNoSolvent") {
		REQUIRE(doEightResiduesNoSolvent(envmode));		
	}
}


namespace StressTesting {
	TEST_CASE("StressTesting::RepeatPool50x") {
		doPool50x(envmode);
	}
}