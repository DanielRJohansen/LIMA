#pragma once

#include "Environment.h"
#include "Printer.h"
#include "Utilities.h"
#include "LimaTypes.cuh"
#include "Filehandling.h"

#include <iostream>
#include <string>
#include <algorithm>
#include <iostream>
#include <functional>
#include <filesystem>

namespace TestUtils {
#ifndef __linux__
	const fs::path simulations_dir = "C:/Users/Daniel/git_repo/LIMA_data/";
#else
	const fs::path simulations_dir = "/home/lima/Downloads/LIMA_data/";
#endif

	fs::path getMostSuitableGroFile(const fs::path& workdir) {
		const fs::path em = workdir / "molecule/em.gro";
		const fs::path conf = workdir / "molecule/conf.gro";
		if (std::filesystem::exists(em)) {
			return em;
		}
		else {
			return conf;
		}
	}

	void CleanDirectory(const fs::path& dir) {
		if (dir.string().find("LIMA_data") == std::string::npos) {
			throw std::runtime_error("LIMA is not allowed to clean this directory");
		}

		if (!fs::exists(dir) || !fs::is_directory(dir)) {
			std::cerr << "Path does not exist or is not a directory: " << dir << std::endl;
			return;
		}

		for (const auto& entry : fs::directory_iterator(dir)) {
			if (fs::is_regular_file(entry.path())) {
				fs::remove(entry.path());  // Delete the file
			}
		}
	}

	void CleanDirIfNotContains(const fs::path& dir, const std::string& except) {
		if (dir.string().find("LIMA_data") == std::string::npos) {
			throw std::runtime_error("LIMA is not allowed to clean this directory");
		}

		if (!fs::exists(dir) || !fs::is_directory(dir)) {
			std::cerr << "Path does not exist or is not a directory: " << dir << std::endl;
			return;
		}

		// Remove all files that do not contain "reference" in their name
		for (auto& p : fs::recursive_directory_iterator(dir)) {
			if (fs::is_regular_file(p) && p.path().filename().string().find(except) == std::string::npos) {
				try {
					fs::remove(p);
				}
				catch (const fs::filesystem_error& e) {
					std::cerr << "Error removing " << p.path() << ": " << e.what() << '\n';
				}
			}
		}

		// Collect all directories in a vector
		std::vector<fs::path> directories;
		for (auto& p : fs::recursive_directory_iterator(dir)) {
			if (fs::is_directory(p)) {
				directories.push_back(p);
			}
		}

		// Sort directories in reverse order to remove from the deepest level first
		std::sort(directories.rbegin(), directories.rend());

		// Remove empty directories
		for (const auto& p : directories) {
			if (fs::is_empty(p)) {
				try {
					fs::remove(p);
				}
				catch (const fs::filesystem_error& e) {
					std::cerr << "Error removing " << p << ": " << e.what() << '\n';
				}
			}
		}
	}

	// Creates a simulation from the folder which should contain a molecule with conf and topol
	// Returns an environment where solvents and compound can still be modified, and nothing (i hope) have
	// yet been moved to device. I should find a way to enforce this...
	static std::unique_ptr<Environment> basicSetup(const std::string& foldername, LAL::optional<SimParams> simparams, EnvMode envmode) {
		
		const fs::path work_folder = simulations_dir / foldername;
		const GroFile conf{getMostSuitableGroFile(work_folder)};
		const TopologyFile topol {work_folder / "molecule/topol.top"};
		const fs::path simpar = work_folder / "sim_params.txt";

		auto env = std::make_unique<Environment>(work_folder, envmode);

		const SimParams ip = simparams.hasValue()
			? simparams.value()
			: SimParams{ simpar };
		

		env->CreateSimulation(conf, topol, ip);

		return std::move(env);
	}

	// assumes that all the values are positive
	bool isOutsideAllowedRange(float value, float target, float maxError=0.1) {
		if (isnan(value)) 
			return true;

		const float error = std::abs(value - target) / target;
		return error > maxError;
	}

	bool isAboveVcThreshold(float value, float target) {
		return value > target;
	}


	bool CompareVecWithFile(const std::vector<Float3>& vec, const fs::path& path, float errorThreshold, bool overwriteFile) {
		if (overwriteFile) {
			FileUtils::WriteVectorToBinaryFile(path, vec);
			return true;
		}

		const std::vector<Float3> fileVec = FileUtils::ReadBinaryFileIntoVector<Float3>(path);
		if (vec.size() != fileVec.size()) {
			return false;
		}
		
		for (size_t i = 0; i < vec.size(); i++) {
			if ((vec[i] - fileVec[i]).len() > errorThreshold) {
				return false;
			}
		}

		return true;
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
		LimaUnittestResult( bool success, const std::string err, const bool print_now) :
			success(success),
			error_description(!err.empty() ? err : success ? "Success" : "Fail")
		{
			if (print_now) {
				printStatus();
			}
		}


		void printStatus() const {
			if (success) {
				setConsoleTextColorGreen();
			}
			else {
				setConsoleTextColorRed();
			}


			if (error_description.length() > 55) { std::cout << "\n\t"; }
			std::cout << error_description << "\n";


			setConsoleTextColorDefault();
		}

		bool success;
		std::string error_description;		
	};

#define ASSERT(condition, errorMsg) \
    do { \
        if (!(condition)) { \
            std::string msg = errorMsg; \
            return LimaUnittestResult{ false, msg, (envmode) == Full }; \
        } \
    } while (0)

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
				testresult = std::make_unique<LimaUnittestResult>(LimaUnittestResult{ false, err_desc, true });
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
			if (successCount == tests.size()) {
				setConsoleTextColorGreen();
			}
			else {
				setConsoleTextColorRed();
			}
			
			std::printf("\n\n#--- Unittesting finished with %d successes of %zu tests ---#\n\n", successCount, tests.size());

			for (const auto& test : tests) {
				if (!test->testresult->success) {
					test->testresult->printStatus();
				}
			}

			setConsoleTextColorDefault();
		}

		void addTest(std::unique_ptr<LimaUnittest> test) {
			test->execute();

			if (test->testresult->success) { successCount++; }

			tests.push_back(std::move(test));
		}

	private:
		std::vector<std::unique_ptr<LimaUnittest>> tests;
		int successCount = 0;
	};




	static LimaUnittestResult loadAndRunBasicSimulation(
		const string& folder_name,
		EnvMode envmode,
		float max_vc = 0.001,
		float max_gradient=1e-7,
		LAL::optional<SimParams> ip = {}
	)
	{		
		auto env = TestUtils::basicSetup(folder_name, ip, envmode);
		env->run();

		const auto analytics = env->getAnalyzedPackage();
		
		float varcoff = analytics.variance_coefficient;
		

		if (envmode != Headless) {
			analytics.Print();
			//LIMA_Print::printPythonVec("potE", std::vector<float>{ analytics.pot_energy});
			//LIMA_Print::printPythonVec("kinE", std::vector<float>{ analytics.kin_energy});
			//LIMA_Print::printPythonVec("totE", std::vector<float>{ analytics.total_energy});
			//LIMA_Print::plotEnergies(analytics.pot_energy, analytics.kin_energy, analytics.total_energy);
		}
		ASSERT(env->getSimPtr()->simsignals_host.critical_error_encountered == false, "Critical error encountered");
		ASSERT(env->getSimPtr()->getStep() == env->getSimPtr()->simparams_host.n_steps, std::format("Simulation did not finish {}/{}",
			env->getSimPtr()->getStep(), env->getSimPtr()->simparams_host.n_steps));

		const auto result = evaluateTest({ varcoff }, max_vc, {analytics.energy_gradient}, max_gradient);

		return LimaUnittestResult{ result.first, result.second, envmode == Full };
	}

	void stressTest(std::function<void()> func, size_t reps) {
		for (size_t i = 0; i < reps; i++) {
			func();
		}
	}

	string compareFilesBitwise(const std::filesystem::path& path1, const std::filesystem::path& path2) {
		// Open the files
		std::ifstream file1(path1, std::ifstream::ate);
		std::ifstream file2(path2, std::ifstream::ate);

		// Check if both files are open
		if (!file1.is_open() || !file2.is_open()) {
			return std::format("Failed to open either or both files \n\t\t{} \n\t\t{}", path1.string(), path2.string());
		}

		// Validate the files. If they are not even 50 bytes long, something is surely wrong
		if (file1.tellg() < 50 || file2.tellg() < 50) {
			return std::format("Expected files to be atleast 50 bytes long \n\t\t{} \n\t\t{}", path1.string(), path2.string());
		}
		//// Compare file sizes
		//file1.seekg(0, std::ifstream::end);
		//file2.seekg(0, std::ifstream::end);
		//if (file1.tellg() != file2.tellg()) {
		//	return "Files are of different length";
		//}
		// 

		// Move ptr back to beginning of file
		file1.seekg(0, std::ifstream::beg);
		file2.seekg(0, std::ifstream::beg);

		// Compare the contents
		if (!std::equal(std::istreambuf_iterator<char>(file1.rdbuf()), std::istreambuf_iterator<char>(), std::istreambuf_iterator<char>(file2.rdbuf()))) {
			return std::format("Files did not match bit for bit \n\t\t{} \n\t\t{}", path1.string(), path2.string());
		};
		return "";
	}


	void CompareForces1To1(const fs::path& workDir, Environment& env, bool overwriteRef) {
		const ParticleDataBuffer<Float3>* forcebuffer = env.getSimPtr()->forceBuffer.get();
		std::vector<Float3> forces(forcebuffer->GetBufferAtStep(0), forcebuffer->GetBufferAtStep(0) + forcebuffer->n_particles_upperbound);

		if (overwriteRef)
			FileUtils::WriteVectorToBinaryFile(workDir/ "forces.bin", forces);

		const std::vector<Float3> forcesRef = FileUtils::ReadBinaryFileIntoVector<Float3>(workDir / "forces.bin");

		std::vector<float> errors(forces.size()); // Pre-allocate the vector
		std::transform(forces.begin(), forces.end(), forcesRef.begin(), errors.begin(),
			[](const Float3& a, const Float3& b) { return (a - b).len(); });



		FileUtils::WriteVectorToBinaryFile(workDir / "errors.bin", errors);
		std::string command = "python " + (FileUtils::GetLimaDir() / "dev/PyTools/pdf.py").string() + " \"" + (workDir / "errors.bin").string() + "\"";
		std::system(command.c_str());
	}

} // namespace TestUtils


