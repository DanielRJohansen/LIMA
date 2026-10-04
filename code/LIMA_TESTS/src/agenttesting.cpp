#include "Tests.h"
#include "Format.h"
#include "Display.h"
#include "DisplayTests.h"
#include "EnergyMinimizationtests.h"
#include "Environment.h"
#include "Programs.h"
#include "TestUtils.h"
#include "Workflow.h"

#include <algorithm>
#include <chrono>
#include <cfloat>
#include <cmath>
#include <filesystem>
#include <format>
#include <iostream>
#include <memory>
#include <string_view>
#include <thread>
#include <vector>


namespace {
	void WriteTestSystem(const std::filesystem::path& directory, GroFile coordinates,
		TopologyFile topology, bool solvate,
		int solventDensity = SimulationBuilder::defaultSolventsPerNm3) {
		std::filesystem::create_directories(directory);
		if (solvate) {
			const float boxLength = std::ceil(std::max({
				coordinates.box_size.x, coordinates.box_size.y, coordinates.box_size.z }));
			coordinates.box_size = Float3{ boxLength };
			std::cout << "Solvating " << directory.filename() << " in a " << boxLength
				<< " nm box at density " << solventDensity << "\n";
			SimulationBuilder::SolvateGrofile(coordinates, topology, solventDensity);
		}
		coordinates.printToFile(directory / "conf.gro");
		topology.printToFile(directory / "topol.top");
	}

	void PlaceInPaddedBox(GroFile& coordinates, float padding) {
		Float3 minimum{ FLT_MAX };
		Float3 maximum{ -FLT_MAX };
		for (const auto& atom : coordinates.atoms) {
			minimum = Float3::ElementwiseMin(minimum, atom.position);
			maximum = Float3::ElementwiseMax(maximum, atom.position);
		}
		const Float3 offset = Float3{ padding } - minimum;
		for (auto& atom : coordinates.atoms)
			atom.position += offset;
		coordinates.box_size = maximum - minimum + Float3{ 2.f * padding };
	}

	void PrepareEnergyMinimizationTests() {
		const std::filesystem::path testRoot = TestUtils::HeavyTestsDir() / "EnergyMinimizationTests";

		const auto MakeMembrane = [&](const std::string& name, Float3 boxSize,
			const MembraneGeometry::Figure& geometry, bool solvate) {
			const std::filesystem::path directory = testRoot / name;
			Lipids::Selection lipids{ Lipids::Select{ "DMPC", directory, 100. } };
			GroFile coordinates;
			coordinates.box_size = boxSize;
			coordinates.title = name;
			TopologyFile topology;
			topology.SetSystem(name);
			SimulationBuilder::CreateMembrane(coordinates, topology, lipids, geometry, 42);
			WriteTestSystem(directory, std::move(coordinates), std::move(topology), solvate);
		};

		MakeMembrane("membrane_plane_unsolvated", Float3{ 8.f }, MembraneGeometry::Plane{ 4.f }, false);
		MakeMembrane("membrane_plane_solvated", Float3{ 8.f }, MembraneGeometry::Plane{ 4.f }, true);
		MakeMembrane("membrane_sphere_unsolvated", Float3{ 16.f },
			MembraneGeometry::Sphere{ Float3{ 8.f }, 5.f }, false);
		MakeMembrane("membrane_sphere_solvated", Float3{ 16.f },
			MembraneGeometry::Sphere{ Float3{ 8.f }, 5.f }, true);
		MakeMembrane("membrane_ellipsoid_unsolvated", Float3{ 20.f },
			MembraneGeometry::Ellipsoid{ Float3{ 10.f }, Float3{ 5.f, 6.f, 7.f } }, false);
		MakeMembrane("membrane_ellipsoid_solvated", Float3{ 20.f },
			MembraneGeometry::Ellipsoid{ Float3{ 10.f }, Float3{ 5.f, 6.f, 7.f } }, true);

		auto t4 = Programs::ToGmx(TestUtils::AutomatedTestsDir() / "pdb2gmx" / "6lzm.pdb");
		WriteTestSystem(testRoot / "t4_solvated", std::move(t4.grofile), std::move(t4.topology), true);

		const std::filesystem::path stmvSource = TestUtils::HeavyTestsDir() / "benchmarking" / "stmv";
		const std::filesystem::path stmvTarget = testRoot / "stmv_solvated";
		std::filesystem::create_directories(stmvTarget);
		std::filesystem::copy_file(stmvSource / "conf.gro", stmvTarget / "conf.gro",
			std::filesystem::copy_options::overwrite_existing);
		std::filesystem::copy_file(stmvSource / "topol.top", stmvTarget / "topol.top",
			std::filesystem::copy_options::overwrite_existing);
		for (const auto& entry : std::filesystem::directory_iterator(stmvSource)) {
			if (entry.is_regular_file() && entry.path().extension() == ".itp")
				std::filesystem::copy_file(entry.path(), stmvTarget / entry.path().filename(),
					std::filesystem::copy_options::overwrite_existing);
		}

		auto conversion = Programs::ToGmx(TestUtils::HeavyTestsDir() / "fileconversions" / "3J3Q.cif");
		PlaceInPaddedBox(conversion.grofile, 1.f);
		size_t moleculeIndex = 0;
		for (auto& [_, molecule] : conversion.topology.moleculetypes)
			molecule->includePath = Lima::Format("molecule_{:04}.itp", moleculeIndex++);
		WriteTestSystem(testRoot / "3j3q_solvated", std::move(conversion.grofile),
			std::move(conversion.topology), true, 1);
		std::cout << "Prepared non-minimized systems in " << testRoot << '\n';
	}

	int RunDisplayPreview(bool tilePreview) {
		const GroFile molecule{ TestUtils::AutomatedTestsDir() / "T4Lysozyme" / "molecule" / "conf.gro" };
		Display display;
		display.allowUserInputs = !tilePreview;
		for (int i = 0; i < (tilePreview ? 9 : 4); ++i) {
			auto task = std::make_unique<Rendering::AtomRenderTask>(molecule, false);
			task->simStatus.step = 24000 + i * 1000;
			task->simStatus.temperature = 300.12f + i;
			task->simStatus.maxForce = 1.23e3f;
			task->simStatus.expectedTimeToFinish = 154.;
			task->simStatus.avgStepTime = .842f;
			task->simStatus.simulationPerformance = 205.23f;
			display.Submit(i, std::move(task), false, nullptr, "Preview " + std::to_string(i + 1));
		}
		while (!display.DisplaySelfTerminated())
			std::this_thread::sleep_for(std::chrono::milliseconds(50));
		if (display.displayThreadException)
			std::rethrow_exception(display.displayThreadException);
		return 0;
	}

	void PrintUsage() {
		std::cout
			<< "Usage: agenttesting OPTION\n"
			<< "  --prepare-energy-minimization-tests\n"
			<< "  --energy-minimization-tests [DIRECTORY]\n"
			<< "  --display-tiles-tests\n"
			<< "  --display-preview\n"
			<< "  --display-tiles-preview\n"
			<< "  --environment-batch-tests\n"
			<< "  --engine-batch-tests\n"
			<< "  --showcase-multisim\n";
	}
}

int main(int argc, char** argv) {
	try {
		if (argc == 2 && std::string_view(argv[1]) == "--prepare-energy-minimization-tests") {
			PrepareEnergyMinimizationTests();
			return 0;
		}
		if (argc >= 2 && std::string_view(argv[1]) == "--energy-minimization-tests") {
			if (argc > 3)
				throw std::runtime_error("--energy-minimization-tests accepts at most one directory");
			const std::filesystem::path testRoot = argc == 3
				? std::filesystem::path{ argv[2] }
				: TestUtils::HeavyTestsDir() / "EnergyMinimizationTests";
			return EnergyMinimizationTests::Run(testRoot);
		}
		if (argc == 2 && std::string_view(argv[1]) == "--display-tiles-tests") {
			DisplayTests::Run();
			return 0;
		}
		if (argc == 2 && std::string_view(argv[1]) == "--display-preview")
			return RunDisplayPreview(false);
		if (argc == 2 && std::string_view(argv[1]) == "--display-tiles-preview")
			return RunDisplayPreview(true);
		if (argc == 2 && std::string_view(argv[1]) == "--environment-batch-tests") {
			BatchingTests::RunSchedulerTests(Environment::Get());
			return 0;
		}
		if (argc == 2 && std::string_view(argv[1]) == "--engine-batch-tests") {
			EngineBatchTests::RunAll();
			return 0;
		}
		if (argc == 2 && std::string_view(argv[1]) == "--showcase-multisim") {
			EngineBatchTests::ShowcaseMultisim();
			return 0;
		}

		PrintUsage();
		return 1;
	}
	catch (const std::exception& ex) {
		std::cerr << ex.what() << '\n';
		return 1;
	}
}
