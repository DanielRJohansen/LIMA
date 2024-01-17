#pragma once

#include "Environment.h"
#include "Printer.h"
#include "Utilities.h"
#include "LimaTypes.cuh"
#include "EngineUtils.cuh"


#include <iostream>
#include <string>
#include <algorithm>

#include <iostream>
#include <optional>
#include <functional>
#include <filesystem>

namespace TestUtils {
#ifndef __linux__
	//const std::string simulations_dir = "C:/PROJECTS/Quantom/Simulation/";


	const std::string simulations_dir = "C:/Users/Daniel/git_repo/LIMA_data/";
#else
	const std::string simulations_dir = "/home/lima/Desktop/LIMA/Simulations/";
#endif

	std::string getMostSuitableGroFile(const std::string& workdir) {
		const std::string em = workdir + "/molecule/em.gro";
		const std::string conf = workdir + "/molecule/conf.gro";
		if (std::filesystem::exists(em)) {
			return em;
		}
		else {
			return conf;
		}
	}

	// Creates a simulation from the folder which should contain a molecule with conf and topol
	// Returns an environment where solvents and compound can still be modified, and nothing (i hope) have
	// yet been moved to device. I should find a way to enforce this...
	static std::unique_ptr<Environment> basicSetup(const std::string& foldername, LAL::optional<SimParams> simparams, EnvMode envmode) {
		
		const std::string work_folder = simulations_dir + foldername;
		//const std::string conf = work_folder + "molecule/conf.gro";
		const std::string conf = getMostSuitableGroFile(work_folder);
		const std::string topol = work_folder + "/molecule/topol.top";
		const std::string simpar = work_folder + "/sim_params.txt";

		auto env = std::make_unique<Environment>(work_folder, envmode, false);

		const SimParams ip = simparams.hasValue()
			? simparams.value()
			: env->loadSimParams(simpar);

		env->CreateSimulation(conf, topol, ip);

		return std::move(env);
	}

	// assumes that all the values are positive
	bool isOutsideAllowedRange(float value, float target, float allowedRangeFromTarget=0.05) {
		if (isnan(value)) 
			return true;

		const float error = std::abs(value - target);
		return error > target * allowedRangeFromTarget;
	}

	bool isAboveVcThreshold(float value, float target) {
		return value > target;
	}


	/// <summary></summary>	
	/// <returns>{success, error_string(empty if successful)}</returns>
	std::pair<bool, std::string> evaluateTest(std::vector<float> VCs, float target_vc, std::vector<float> energy_gradients, float max_energygradient_abs)
	{
		// Pick the correct evaluate function depending on if we have multiple VCs. Cant set a target vc to keep, if we have different sims ;)
		auto evaluateVC = [&](float vc) {
			if (VCs.size() > 1) {
				return isAboveVcThreshold(vc, target_vc);
			}
			else {
				return isOutsideAllowedRange(vc, target_vc);
			}
		};


		for (auto& vc : VCs) {
			if (evaluateVC(vc)) {
				return { false, std::format("Variance Coefficient of {:.3e} was too far from the target {:.3e}", vc, target_vc) };
			}
		}

		for (auto& gradient : energy_gradients) {
			if (isnan(gradient) || abs(gradient) > max_energygradient_abs) {
				return { false, std::format("Energygradient of {:.3e} superceeded the max of {:.3e}", gradient, max_energygradient_abs) };
			}
		}

		float highest_vc = *std::max_element(VCs.begin(), VCs.end());
		return { true, std::format("VC {:.3e} / {:.3e}", highest_vc, target_vc)};
	}

	static void setConsoleTextColorRed() { std::cout << "\033[31m"; }
	static void setConsoleTextColorGreen() { std::cout << "\033[32m"; }
	static void setConsoleTextColorDefault() { std::cout << "\033[0m"; }

	struct LimaUnittestResult {
		enum TestStatus { SUCCESS, FAIL, THROW };

		LimaUnittestResult( TestStatus status, const std::string err, const bool print_now) :
			status(status),
			error_description(err)
		{
			if (print_now) {
				printStatus();
			}
		}


		void printStatus() const {
			std::string status_str;// = status == SUCCESS ? "Success" : "Failure";

			if (status == SUCCESS) {
				status_str = "Success";
				setConsoleTextColorGreen();
			}
			else {
				status_str = "Fail"   ;
				setConsoleTextColorRed();
			}

			std::cout << " status: " << status_str;

			if (error_description.length() > 30) { std::cout << "\n"; }
			std::cout << "\t" << error_description << "\n";


			setConsoleTextColorDefault();
		}

		bool success() const { return status == SUCCESS; }

		TestStatus status;		
		std::string error_description;		
	};

	struct LimaUnittest {
		LimaUnittest(const std::string& name, std::function<LimaUnittestResult()> test) :
			name(name),
			test(test)
		{}

		void execute() {
			try {
				std::cout << "Test " << name << " ";
				testresult = std::make_unique<LimaUnittestResult>(test());

				int str_len = 6 + name.length();
				while (str_len++ < 61) { std::cout << " "; }

				testresult->printStatus();
			}
			catch (const std::runtime_error& ex) {
				const std::string err_desc = "Test threw exception: " + std::string(ex.what());
				testresult = std::make_unique<LimaUnittestResult>(LimaUnittestResult{ LimaUnittestResult::THROW, err_desc, true });
			}
		}

		const std::function<LimaUnittestResult()> test;
		std::unique_ptr<LimaUnittestResult> testresult;
		const std::string name;
	};


	class LimaUnittestManager {
	public:
		LimaUnittestManager(){}
		~LimaUnittestManager() {
			if (successes == tests.size()) {
				setConsoleTextColorGreen();
			}
			else {
				setConsoleTextColorRed();
			}
			
			std::printf("\n\n#--- Unittesting finished with %d successes of %zu tests ---#\n\n", successes, tests.size());

			for (const auto& test : tests) {
				if (!test->testresult->success()) {
					test->testresult->printStatus();
				}
			}

			setConsoleTextColorDefault();
		}

		void addTest(std::unique_ptr<LimaUnittest> test) {
			test->execute();

			if (test->testresult->success()) { successes++; }

			tests.push_back(std::move(test));
		}

	private:
		std::vector<std::unique_ptr<LimaUnittest>> tests;
		int successes = 0;
	};




	static LimaUnittestResult loadAndRunBasicSimulation(
		const string& folder_name,
		EnvMode envmode,
		float max_vc = 0.001,
		float max_gradient=1e-7,
		LAL::optional<SimParams> ip = {},
		bool em_variant = false
	)
	{
		auto env = TestUtils::basicSetup(folder_name, ip, envmode);
		env->run(em_variant);

		const auto analytics = env->getAnalyzedPackage();
		
		float varcoff = analytics->variance_coefficient;
		

		if (envmode != Headless) {
			Analyzer::printEnergy(analytics);
			LIMA_Print::printMatlabVec("varcoffs", std::vector<float>{ varcoff });
			LIMA_Print::printMatlabVec("energy_gradients", std::vector<float>{ analytics->energy_gradient });
		}


		const auto result = evaluateTest({ varcoff }, max_vc, {analytics->energy_gradient}, max_gradient);
		const auto status = result.first == true ? LimaUnittestResult::SUCCESS : LimaUnittestResult::FAIL;

		return LimaUnittestResult{ status, result.second, envmode == Full };
	}

	void stressTest(std::function<void()> func, size_t reps) {
		for (size_t i = 0; i < reps; i++) {
			func();
		}
	}
}


