#pragma once

// Declarations of all tests. Each test file is its own .cpp, so editing one test (or main.cpp, to pick which
// test to run) only recompiles that file

#include "TestUtils.h"
#include "Programs.h"
#include "MDFiles.h"
#include "Bodies.cuh"

#include <filesystem>
#include <optional>
#include <string>
#include <vector>

using TestUtils::TestRoutine;
using TestUtils::LimaUnittestResult;

// AlgorithmTests.cpp
namespace KernelAlgorithms {
	LimaUnittestResult WarpSort64_Unittest(EnvMode envmode);
}

// BatchingTests.cpp
namespace BatchingTests {
	SimulationJob MakeT4Job(EnvMode envmode, bool mustRunAlone);
	bool HasIdenticalCoordinates(const Simulation& lhs, const Simulation& rhs);
	TestRoutine TestFourT4BatchMatchReference(Environment& environment, EnvMode envmode);
	SimulationJob MakeSmallJob(int particles, int steps, float dt, bool electrostatics);
	std::vector<SimulationHandle> SubmitTogether(Environment& environment, std::vector<SimulationJob> jobs);
	void RunSchedulerTests(Environment& environment);
}

// Benchmarks.cpp
namespace Benchmarks {
	namespace fs = std::filesystem;
	constexpr int automatedTestRuns = 3;
	template<typename Duration>
	struct PerformanceBounds {
		Duration min;
		Duration max;
	};
	const fs::path TestsDir();
	TestRoutine ToGmxLargeCif(Environment&, EnvMode envmode);
	TestRoutine Load3J3Q(Environment& environment, EnvMode envmode);
	TestRoutine LoadT4(Environment& environment, EnvMode envmode);
	TestRoutine Bench(Environment& environment, EnvMode envmode, fs::path workDir,
			fs::path groPath, fs::path topPath, fs::path simParamsPath,
			PerformanceBounds<std::chrono::microseconds> allowedTimePerStep, int nSteps, int nRuns, int warmupSteps = 0);
	TestRoutine STMV(Environment& environment, EnvMode envmode, int nSteps, int nRuns = 1);
	TestRoutine T4(Environment& environment, EnvMode envmode, int nSteps = 500, int nRuns = 1);
}

// ElectrostaticsTests.cpp
namespace ElectrostaticsTests {
	TestRoutine CoulombForceSanityCheck(Environment&, EnvMode envmode);
	Float3 GetPositionOfParticleRelativeToSelfUsingTheWierdLogicOfTheKernel(const Float3& posOtherAbs, const NodeIndex nodeindexSelf);
	LimaUnittestResult TestAttractiveParticlesInteractingWithESandLJ(EnvMode envmode);
	void MakeChargeParticlesSim(
			GroFile& grofile, TopologyFile& topfile, const fs::path& workDir,
			const float boxLen, const AtomsSelection& atomsSelection, float particlesPerNm3);
	TestRoutine TestChargedParticlesVelocityInUniformElectricField(
			Environment& environment, EnvMode envmode);
	TestRoutine TestLongrangeEsNoLJTwoParticles(
			Environment& environment, EnvMode envmode);
	LimaUnittestResult PlotPmePotAsFactorOfDistance(EnvMode envmode);
	LimaUnittestResult TestConsistentEnergyWhenGoingFromLresToSres(EnvMode envmode);
	TestRoutine TestLongrangeEsNoLJManyParticles(
			Environment& environment, EnvMode envmode);
}

// EngineBatchTests.cpp
namespace EngineBatchTests {
	void Require(bool value, const char* message);
	std::unique_ptr<Simulation> MakeSimulation(int particles, int steps, float dt, bool electrostatics, bool em = false, int loggingInterval = 2);
	void Run(Simulation& sim);
	void Compare(Simulation& a, Simulation& b);
	void RunAll();
	void ShowcaseMultisim();
}

// FileTests.cpp
namespace FileTests {
	namespace fs = std::filesystem;
	LimaUnittestResult TestFilesAreCachedAsBinaries(EnvMode envmode);
	TestRoutine TestTopologyPreprocessor(Environment&, EnvMode envmode);
}

// ForceComparisons.cpp
namespace ForceComparisons {
	TestRoutine DoAllForceComparisons(Environment& environment, EnvMode envmode);
}

// ForceCorrectness.cpp
namespace ForceCorrectness {
	const fs::path TestsDir();
	SimulationJob MakeJob(const fs::path& workDir, EnvMode envmode, bool analyze = false);
	TestRoutine FinishStabilityTest(Environment& environment, std::string name, EnvMode envmode,
			SimulationJob job);
	TestRoutine doPoolBenchmark(Environment& environment, EnvMode envmode);
	TestRoutine SinglebondForceAndPotentialSanityCheck(Environment& environment, EnvMode envmode);
	TestRoutine SinglebondOscillationTest(Environment& environment, EnvMode envmode);
	struct ExpectedForceEnergy { Float3 force{}; float potential = 0.f; };
	TestRoutine UreyBradleyForceAndPotentialSanityCheck(Environment& environment, EnvMode envmode);
	TestRoutine PairbondForceAndPotentialSanityCheck(Environment& environment, EnvMode envmode);
	TestRoutine doPoolCompSolBenchmark(Environment& environment, EnvMode envmode);
	TestRoutine doSinglebondBenchmark(Environment& environment, EnvMode envmode);
	TestRoutine doAnglebondBenchmark(Environment& environment, EnvMode envmode);
	TestRoutine doDihedralbondBenchmark(Environment& environment, EnvMode envmode);
	TestRoutine doImproperDihedralBenchmark(Environment& environment, EnvMode envmode);
}
namespace VerletintegrationTesting {
	TestRoutine TestIntegration(Environment& environment, EnvMode envmode);
}

// ForcefieldTests.cpp
void ParseForcefieldFromTpr(const std::string& filePath, std::vector<SingleBond::Parameters>& bonds, std::vector<AngleUreyBradleyBond::Parameters>& angles,
    std::vector<DihedralBond::Parameters>& dihedrals, std::vector<ImproperDihedralBond::Parameters>& improperDihedrals, std::vector<std::array<int, 4>>& dihIds);
void ParseForcefieldFromItp(
    const std::string& filePath, 
    std::vector<SingleBond::Parameters>& bonds, 
    std::vector<AngleUreyBradleyBond::Parameters>& angles,
    std::vector<DihedralBond::Parameters>& dihedrals, 
    std::vector<ImproperDihedralBond::Parameters>& improperDihedrals,
    std::vector<std::array<std::string, 2>>& atomnamesSinglebonds,
    std::vector<std::array<std::string, 3>>& atomnamesAnglebonds,
    std::vector<std::array<std::string, 4>>& atomnamesDihedralbonds,
    std::vector<std::array<std::string, 4>>& atomnamesImproperDihedralbonds);
TestRoutine TestLimaChosesSameBondparametersAsGromacs(Environment&, EnvMode envmode);

// MDStability.cpp
namespace TestMDStability {
	SimulationJob MakeEnergyMinJob(const fs::path& workDir, EnvMode envmode);
	TestRoutine LoadEnergyMinAndRunBasicSimulation(
			Environment& environment, EnvMode envmode, std::string folderName, std::string testName);
	TestRoutine TestDeterministic(Environment& environment, EnvMode envmode);
	bool doMoleculeTranslationTest(std::string foldername);
}

// MembraneBuilder.cpp
namespace TestMembraneBuilder {
	namespace fs = std::filesystem;
	std::vector<Float3> ResidueCenters(const GroFile& grofile);
	float NearestNeighborSpacingVariation(const std::vector<Float3>& directions);
	float MinimumNearestNeighborSpacing(const std::vector<Float3>& points);
	TestRoutine TestSphericalMembraneBuilder(Environment&, EnvMode envmode);
	TestRoutine TestBuildmembraneSmall(Environment& environment, EnvMode envmode, bool do_em);
	TestRoutine TestBuildmembraneWithCustomlipidAndCustomForcefield(Environment& environment, EnvMode envmode);
	TestRoutine TestAllStockholmlipids(Environment& environment, EnvMode envmode);
	LimaUnittestResult BuildAndRelaxVesicle(EnvMode envmode);
}

// MinorPrograms.cpp
namespace TestMinorPrograms {
	namespace fs = std::filesystem;
	LimaUnittestResult InsertMoleculesAndDoStaticbodyEM(EnvMode envmode);
}

// ProgramsTests.cpp
namespace ProgramsTests {
	LimaUnittestResult TestBuildMembrane(EnvMode envmode);
	TestRoutine TestToGmx_pdbfile(Environment&, EnvMode envmode);
	TestRoutine TestToGmx_ciffile(Environment&, EnvMode envmode);
	TestRoutine TestToGmx_multichain(Environment&, EnvMode envmode);
}

// SetupTests.cpp
TestRoutine TestBoxIsSavedCorrectlyBetweenSimulations(Environment& environment, EnvMode envmode);
