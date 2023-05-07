#include "LIMA_TESTS/include/ForceCorrectness.h"
#include "LIMA_TESTS/include/MDStability.h"



//
//// Due to the following define, i belive all includes to test-headers must come after
//#define CATCH_CONFIG_MAIN
//#include "Catch2/catch.hpp"	// Don't include windows.h after this line
//
//#include "LIMA_TESTS/include/MDStability.h"
//
//
//
//
//int doThing(int i) {
//	printf("Hello world");
//	return i * 2;
//}
//
//TEST_CASE("Internal", "hey") {
//	REQUIRE(doThing(4) == 8);
//}
//
////bool coordPrecesionBenchmark() {
////	Float3 pos1{ 3.5, 4, 4 };
////	Float3 pos2{ 4, 4, 4 };
////
////	Coord c1{ pos1 };
////	Coord c2{ pos2 };
////
////	double dist_c = sqrt(c1.distSqAbs(&c2));
////	double dist_f = (pos1 - pos2).len();
////
////	printf("dist_c %.12f dist_f %.12f\n", dist_c, dist_f);
////
////	return true;
////}


int main() {

	/*coordPrecesionBenchmark();
	return 0;*/



	//doPoolBenchmark();

	//basicBenchmark(env);

	//doProteinBenchmark(env);
	//doPoolBenchmark(env, "Pool");			// Two 1-particle molecules colliding
	//doPoolCompSolBenchmark(env, "PoolCompSol");	// One 1-particle molecule colliding with 1 solvent
	//doSpringBenchmark(env);
	//doAngleBenchmark(env);	// Doesn't work currently
	

	//loadAndRunBasicSimulation(env, "TorsionBenchmark");
	//loadAndRunBasicSimulation("Met");
	//loadAndRunBasicSimulation("T4LysozymeNoSolventSmall");
	//loadAndRunBasicSimulation(env, "T4LysozymeNoSolvent");
	testNearestSolventSolventAfterEnergyMinimizationIsDecreasing();
	//loadAndRunBasicSimulation("SolventBenchmark", 0.0002);
	//loadAndRunBasicSimulation(env, "T4Lysozyme");
	//loadAndRunBasicSimulation(env, "4ake");// TOO big, almost 20 nm long!
	
	//doMoleculeTranslationTest("T4LysozymeNoSolventSmall");
	//doMoleculeTranslationTest("T4LysozymeNoSolvent");


	return 0;
}

