#include "EnergyMinimizationtests.h"
#include "Format.h"

#include "Environment.h"
#include "MDFiles.h"
#include "SimParams.h"
#include "Simulation.cuh"

#include <algorithm>
#include <cmath>
#include <format>
#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>
#include <vector>

namespace fs = std::filesystem;

namespace EnergyMinimizationTests {
	namespace {
		constexpr float roughForceTolerance = 1000.f;
		constexpr float fineForceTolerance = 200.f;
		constexpr int maximumSteps = 2000;
		constexpr std::chrono::duration<double> maximumRunTime{ 20. };

		struct TestCase {
			std::string name;
			fs::path directory;
			fs::path coordinates;
			fs::path topology;
		};

		struct Result {
			std::string name;
			std::string error;
			std::optional<int64_t> roughStep;
			std::optional<int64_t> fineStep;
			int64_t finalStep = 0;
			float initialForce = NAN;
			float minimumForce = NAN;
			float finalForce = NAN;
			double engineSeconds = 0.;
			size_t forceIncreaseCount = 0;
			bool finite = true;
		};

		// Optional environment configuration for comparing EM algorithms:
		// LIMA_EM_LABEL names the results files, LIMA_EM_CASES is a comma-separated list of case names to run
		std::string ResultsLabel() {
			const char* label = std::getenv("LIMA_EM_LABEL");
			return label ? std::string{ "_" } + label : std::string{};
		}

		bool IsCaseSelected(const std::string& name) {
			const char* cases = std::getenv("LIMA_EM_CASES");
			if (!cases) return true;
			std::stringstream stream{ cases };
			std::string selected;
			while (std::getline(stream, selected, ','))
				if (selected == name) return true;
			return false;
		}

		void WriteCurve(const fs::path& path, const Simulation& simulation) {
			fs::create_directories(path.parent_path());
			std::ofstream output{ path, std::ios::trunc };
			output << "step,max_force,dt\n";
			for (const auto& entry : simulation.emLog)
				output << Lima::Format("{},{},{}\n", entry.step, entry.maxForce, entry.dt);
		}

		std::vector<TestCase> FindTestCases(const fs::path& testRoot) {
			std::vector<TestCase> testCases;
			if (!fs::is_directory(testRoot))
				throw std::runtime_error("Energy-minimization test directory does not exist: " + testRoot.string());

			for (const auto& entry : fs::directory_iterator(testRoot)) {
				if (!entry.is_directory())
					continue;
				const fs::path coordinates = entry.path() / "conf.gro";
				const fs::path topology = entry.path() / "topol.top";
				if (entry.path().filename() == "3j3q_solvated")
					continue; // Kept as a heavyweight corpus case, but excluded from the quick suite for now.
				if (!IsCaseSelected(entry.path().filename().string()))
					continue;
				if (fs::exists(coordinates) && fs::exists(topology))
					testCases.push_back({ entry.path().filename().string(), entry.path(), coordinates, topology });
			}
			std::ranges::sort(testCases, {}, &TestCase::name);
			return testCases;
		}

		Result RunTestCase(const TestCase& testCase) {
			SimParams params = SimParams::BasicEMSimParams(fineForceTolerance);
			params.n_steps = maximumSteps;
			params.data_logging_interval = 0; // Trajectory buffers are preallocated for n_steps, which does not fit in memory for 3J3Q

			SimulationJob job{
				testCase.directory,
				GroFile{ testCase.coordinates },
				TopologyFile{ testCase.topology },
				params,
				EnvMode::Headless
			};
			job.name = testCase.name;
			job.mustRunAlone = true;
			job.maxRunTime = maximumRunTime;
			SimulationResult simulationResult = Environment::Get().Submit(std::move(job)).Get();

			Result result;
			result.name = testCase.name;
			result.engineSeconds = simulationResult.engineTime.count();
			WriteCurve(testCase.directory.parent_path() / ("curves" + ResultsLabel()) / (testCase.name + ".csv"), *simulationResult.simulation);
			const auto& forces = simulationResult.simulation->maxForceBuffer;
			if (forces.empty()) {
				result.finite = false;
				return result;
			}

			result.initialForce = forces.front().second;
			result.finalForce = forces.back().second;
			result.finalStep = forces.back().first;
			result.minimumForce = result.initialForce;
			for (size_t i = 0; i < forces.size(); ++i) {
				const auto [step, force] = forces[i];
				result.finite = result.finite && std::isfinite(force);
				result.minimumForce = std::min(result.minimumForce, force);
				if (!result.roughStep && force <= roughForceTolerance)
					result.roughStep = step;
				if (!result.fineStep && force <= fineForceTolerance)
					result.fineStep = step;
				if (i > 0 && force > forces[i - 1].second)
					++result.forceIncreaseCount;
			}
			if (simulationResult.execution.timedOut)
				result.error = "time limit reached";
			else if (!result.fineStep)
				result.error = "step limit reached";
			return result;
		}

		std::string StepString(const std::optional<int64_t> step) {
			return step ? std::to_string(*step) : "not reached";
		}

		void WriteResults(const fs::path& path, const std::vector<Result>& results) {
			std::ofstream output{ path, std::ios::trunc };
			if (!output)
				throw std::runtime_error("Could not write energy-minimization results: " + path.string());
			output << "simulation,rough_step_1000,fine_step_200,final_step,initial_max_force,minimum_max_force,final_max_force,engine_seconds,force_increase_count,finite,error\n";
			for (const Result& result : results) {
				output << result.name << ',';
				if (result.roughStep) output << *result.roughStep;
				output << ',';
				if (result.fineStep) output << *result.fineStep;
				output << ',' << result.finalStep << ',' << result.initialForce << ',' << result.minimumForce << ',' << result.finalForce
					<< ',' << result.engineSeconds << ',' << result.forceIncreaseCount << ',' << result.finite << ',' << result.error << '\n';
			}
		}

		std::string EscapeHtml(std::string_view text) {
			std::string escaped;
			for (const char character : text) {
				switch (character) {
				case '&': escaped += "&amp;"; break;
				case '<': escaped += "&lt;"; break;
				case '>': escaped += "&gt;"; break;
				case '"': escaped += "&quot;"; break;
				default: escaped += character; break;
				}
			}
			return escaped;
		}

		std::string ForceString(float force) {
			return std::isfinite(force) ? Lima::Format("{:.2f}", force) : "—";
		}

		void WriteResultsHtml(const fs::path& path, const std::vector<Result>& results) {
			std::ofstream output{ path, std::ios::trunc };
			if (!output)
				throw std::runtime_error("Could not write energy-minimization HTML results: " + path.string());

			output << R"(<!doctype html>
<html lang="en"><head><meta charset="utf-8"><title>Energy minimization results</title>
<style>
body { font: 15px/1.45 system-ui, sans-serif; margin: 2rem; color: #202124; background: #fafafa; }
h1 { margin-bottom: .2rem; } p { color: #5f6368; }
table { border-collapse: collapse; width: 100%; background: white; box-shadow: 0 1px 3px #0002; }
th, td { padding: .65rem .75rem; text-align: right; border-bottom: 1px solid #e6e6e6; } th:first-child, td:first-child, td:last-child { text-align: left; }
th { background: #f1f3f4; white-space: nowrap; } tr.ok { border-left: 5px solid #188038; } tr.timeout { border-left: 5px solid #f9ab00; } tr.failure { border-left: 5px solid #d93025; }
.status { font-weight: 650; } .ok .status { color: #188038; } .timeout .status { color: #b06000; } .failure .status { color: #d93025; }
</style></head><body><h1>Energy minimization baseline</h1>
<p>Targets: rough ≤ 1000 kJ/mol/nm; fine ≤ 200 kJ/mol/nm. Per-case engine time limit: 10 seconds. 3J3Q is excluded from this quick suite.</p>
<table><thead><tr><th>Simulation</th><th>Status</th><th>Steps to 1000</th><th>Steps to 200</th><th>Final step</th><th>Initial max force</th><th>Minimum max force</th><th>Final max force</th><th>Engine time</th><th>Force increases</th><th>Detail</th></tr></thead><tbody>)";
			for (const Result& result : results) {
				const bool converged = result.error.empty() && result.finite && result.fineStep.has_value();
				const bool timedOut = result.error == "time limit reached";
				const std::string rowClass = converged ? "ok" : timedOut ? "timeout" : "failure";
				const std::string status = converged ? "Converged" : timedOut ? "Timed out" : "Failed";
				output << Lima::Format("<tr class=\"{}\"><td>{}</td><td class=\"status\">{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{:.3f} s</td><td>{}</td><td>{}</td></tr>\n",
					rowClass, EscapeHtml(result.name), status, StepString(result.roughStep), StepString(result.fineStep),
					result.finalStep, ForceString(result.initialForce), ForceString(result.minimumForce), ForceString(result.finalForce),
					result.engineSeconds, result.forceIncreaseCount, EscapeHtml(result.error));
			}
			output << "</tbody></table></body></html>\n";
		}
	}

	int Run(const fs::path& testRoot) {
		const std::vector<TestCase> testCases = FindTestCases(testRoot);
		if (testCases.empty())
			throw std::runtime_error("No energy-minimization test cases found in " + testRoot.string());

		std::vector<Result> results;
		results.reserve(testCases.size());
		const fs::path resultPath = testRoot / ("results" + ResultsLabel() + ".csv");
		const fs::path htmlPath = testRoot / ("results" + ResultsLabel() + ".html");
		bool success = true;
		for (const TestCase& testCase : testCases) {
			std::cout << "Energy minimizing " << testCase.name << "...\n";
			Result result;
			try {
				result = RunTestCase(testCase);
				std::cout << Lima::Format("  rough: {}, fine: {}, minimum force: {:.2f}, engine: {:.3f} s\n",
					StepString(result.roughStep), StepString(result.fineStep), result.minimumForce, result.engineSeconds);
			}
			catch (const std::exception& exception) {
				result.name = testCase.name;
				result.finite = false;
				result.error = exception.what();
				std::ranges::replace(result.error, ',', ';');
				std::cout << "  failed: " << result.error << '\n';
			}
			success = success && result.error.empty() && result.finite
				&& result.roughStep.has_value() && result.fineStep.has_value();
			results.push_back(std::move(result));
			WriteResults(resultPath, results);
			WriteResultsHtml(htmlPath, results);
		}

		std::cout << "Wrote " << resultPath << " and " << htmlPath << '\n';
		return success ? 0 : 1;
	}
}
