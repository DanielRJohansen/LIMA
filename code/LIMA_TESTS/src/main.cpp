

#include "LIMA_TESTS/src/ForceCorrectness.h"
#include "LIMA_TESTS/src/MDStability.cuh"
#include "LIMA_TESTS/src/EnergyMinimization.cuh"

constexpr auto envmode = Environment::Mode::Full;
//
//// Due to the following define, i belive all includes to test-headers must come after
//
//
//#include "LIMA_TESTS/include/MDStability.h"
//
//
//
//

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

using namespace TestMDStability;

int main() {

	/*coordPrecesionBenchmark();
	return 0;*/



	//doPoolBenchmark();

	//basicBenchmark(env);

	//doProteinBenchmark(env);
	//doPoolBenchmark(envmode);			// Two 1-particle molecules colliding
	//doPoolCompSolBenchmark(envmode);	// One 1-particle molecule colliding with 1 solvent
	//doSinglebondBenchmark(envmode);
	//doAnglebondBenchmark(envmode);	// Doesn't work currently
	//doDihedralbondBenchmark(envmode);

	//loadAndRunBasicSimulation(env, "TorsionBenchmark");
	//loadAndRunBasicSimulation("Met");
	//loadAndRunBasicSimulation("T4LysozymeNoSolventSmall");
	//loadAndRunBasicSimulation(env, "T4LysozymeNoSolvent");
	//testNearestSolventSolventAfterEnergyMinimizationIsDecreasing();
	loadAndEMAndRunBasicSimulation("SolventBenchmark", envmode, 0.0002);
	//loadAndRunBasicSimulation("T4Lysozyme");
	//loadAndRunBasicSimulation(env, "4ake");// TOO big, almost 20 nm long!
	
	//doMoleculeTranslationTest("T4LysozymeNoSolventSmall");
	//doMoleculeTranslationTest("T4LysozymeNoSolvent");


	return 0;
}

