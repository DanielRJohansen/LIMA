// Due to the following define, i belive all includes to test-headers must come after
#define CATCH_CONFIG_MAIN

#include "LIMA_TESTS/src/ForceCorrectness.h"
#include "LIMA_TESTS/src/MDStability.cuh"

#include <Catch2/catch.hpp>	// Don't include windows.h after this line

constexpr auto envmode = EnvMode::Headless;

namespace TestForceCorrectness {
	using namespace ForceCorrectness;
	TEST_CASE("TestForceCorrectness::Pool") {
		REQUIRE(doPoolBenchmark(envmode).success());			// Two 1-particle molecules colliding
	}

	TEST_CASE("TestForceCorrectness::PoolCarbonSol") {
		REQUIRE(doPoolCompSolBenchmark(envmode).success());
	}

	TEST_CASE("TestForceCorrectness::SingleBond") {
		REQUIRE(doSinglebondBenchmark(envmode).success());
	}

	TEST_CASE("TestForceCorrectness::AngleBond") {
		REQUIRE(doAnglebondBenchmark(envmode).success());
	}

	TEST_CASE("TestForceCorrectness::DihedralBond") {
		REQUIRE(doDihedralbondBenchmark(envmode).success());
	}

	TEST_CASE("TestForceCorrectness::ImproperDihedralBond") {
		REQUIRE(doImproperDihedralBenchmark(envmode).success());
	}

	TEST_CASE("TestForceCorrectness::Methionine") {
		REQUIRE(doMethionineBenchmark(envmode).success());
	}

	TEST_CASE("TestForceCorrectness::TenSolvents") {
		REQUIRE(TestUtils::loadAndRunBasicSimulation("TenSolvents", envmode, 0.0006f).success());
	}
}

namespace TestMDStability {

	TEST_CASE("TestMDStability::EMRunSolvent") {
		REQUIRE(loadAndEMAndRunBasicSimulation("SolventBenchmark", envmode, 0.0002f).success());
	}

	TEST_CASE("TestMDStability::EightResiduesNoSolvent") {
		REQUIRE(doEightResiduesNoSolvent(envmode).success());
	}
}


namespace StressTesting {
	TEST_CASE("StressTesting::RepeatPool50x") {
		doPool50x(envmode);
	}
}