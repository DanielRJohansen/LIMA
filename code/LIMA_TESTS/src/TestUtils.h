#pragma once

#include "LIMA_MD/include/Environment.h"
#include "LIMA_BASE/include/Printer.h"
#include "LIMA_BASE/include/Utilities.h"
#include "LIMA_BASE/include/LimaTypes.cuh"
#include "LIMA_ENGINE/include/EngineUtils.cuh"


#include <iostream>
#include <string>
#include <algorithm>

#include <iostream>
#include <optional>
#include <functional>


namespace TestUtils {

	// Creates a simulation from the folder which should contain a molecule with conf and topol
	// Returns an environment where solvents and compound can still be modified, and nothing (i hope) have
	// yet been moved to device. I should find a way to enforce this...
	static std::unique_ptr<Environment> basicSetup(const std::string& foldername, LAL::optional<InputSimParams> simparams, EnvMode envmode) {
		
		const std::string work_folder = "C:/PROJECTS/Quantom/Simulation/" + foldername + "/";
		const std::string conf = work_folder + "molecule/conf.gro";
		const std::string topol = work_folder + "molecule/topol.top";
		const std::string simpar = work_folder + "sim_params.txt";

		auto env = std::make_unique<Environment>(work_folder, envmode );

		const InputSimParams ip = simparams
			? simparams.value()
			: env->loadInputSimParams(simpar);

		env->CreateSimulation(conf, topol, ip);

		return std::move(env);
	}


	/// <summary></summary>	
	/// <returns>{success, error_string(empty if successful)}</returns>
	std::pair<bool, std::string> evaluateTest(std::vector<float> VCs, float max_vc, std::vector<float> energy_gradients, float max_energygradient_abs) {
		for (auto& vc : VCs) {
			if (isnan(vc) || vc > max_vc) {					
				return { false, std::format("Variance Coefficient of {:.3e} superceeded the max of {:.3e}", vc, max_vc) };
			}
		}

		for (auto& gradient : energy_gradients) {
			if (isnan(gradient) || abs(gradient) > max_energygradient_abs) {
				return { false, std::format("Energygradient of {:.3e} superceeded the max of {:.3e}", gradient, max_energygradient_abs) };
			}
		}

		float highest_vc = *std::max_element(VCs.begin(), VCs.end());
		return { true, std::format("VC {:.3e} / {:.3e}", highest_vc, max_vc)};
	}

	static void setConsoleTextColorRed() { std::cout << "\033[31m"; }
	static void setConsoleTextColorGreen() { std::cout << "\033[32m"; }
	static void setConsoleTextColorDefault() { std::cout << "\033[0m"; }

	struct LimaUnittest {
		enum TestStatus { SUCCESS, FAIL };

		LimaUnittest(const std::string& name, TestStatus status, const std::string err="", bool auto_print = true) :
			name(name), 
			status(status),
			error_description(err),
			auto_print(auto_print) 
		{}
		~LimaUnittest() {
			if (auto_print) { print(); }
		}

		void print() const {
			std::string status_str;// = status == SUCCESS ? "Success" : "Failure";

			if (status == SUCCESS) {
				status_str = "Success";
				setConsoleTextColorGreen();
			}
			else {
				status_str = "Fail";
				setConsoleTextColorRed();
			}

			std::cout << "Test " << name << " status: " << status_str;

			if (error_description.length() > 30) { std::cout << "\n"; }
			std::cout << "\t" << error_description << "\n";


			setConsoleTextColorDefault();
		}

		bool success() const { return status == SUCCESS; }

		const std::string name;				
		const bool auto_print;				// Should the test automatically print it's results
		const TestStatus status;		
		const std::string error_description;		
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
			
			std::printf("\n\n#--- Unittesting finished with %d successes of %d tests ---#\n\n", successes, tests.size());

			for (const auto& test : tests) {
				if (test.status == LimaUnittest::FAIL) {
					test.print();
				}
			}

			setConsoleTextColorDefault();
		}

		void addTest(const LimaUnittest& test) {
			test.print();

			if (test.status == LimaUnittest::SUCCESS) { successes++; }

			tests.push_back(test);
		}

		//void addTest(std::function<LimaUnittest()> test) {
		//	test.print();

		//	if (test.status == LimaUnittest::SUCCESS) { successes++; }

		//	tests.push_back(test);
		//}

	private:
		std::vector<LimaUnittest> tests;
		int successes = 0;
	};




	static LimaUnittest loadAndRunBasicSimulation(
		const string& folder_name,
		EnvMode envmode,
		float max_vc = 0.001,
		float max_gradient=1e-7,
		LAL::optional<InputSimParams> ip = {},
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
		//Analyzer::findAndDumpPiecewiseEnergies(*env->getSimPtr(), env->getWorkdir());


		const auto result = evaluateTest({ varcoff }, max_vc, {analytics->energy_gradient}, max_gradient);
		const auto status = result.first == true ? LimaUnittest::SUCCESS : LimaUnittest::FAIL;

		return LimaUnittest{ "loadAndRunBasicSimulation:" + folder_name, status, result.second, envmode == Full };
	}

	void stressTest(std::function<void()> func, size_t reps) {
		for (size_t i = 0; i < reps; i++) {
			func();
		}
	}
}


