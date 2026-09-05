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
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>

namespace TestUtils {

	fs::path AutomatedTestsDir() { return FileUtils::GetLimaDir() / "tests" / "automatedtests"; }
	fs::path HeavyTestsDir() { return FileUtils::GetLimaDir().parent_path() / "LIMA_data"; }

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

	bool MayModifyDir(const fs::path& path) {
		if (path.string().find("LIMA_data") == std::string::npos && path.string().find("automatedtests") == std::string::npos) {
			throw std::runtime_error("LIMA is not allowed to clean this directory");
			return false;
		}
		return true;
	}

	void TryDeleteFile(const fs::path& path) {
		MayModifyDir(path);
		try {
			fs::remove(path);
		}
		catch (const fs::filesystem_error& e) {
			std::cerr << "Error removing " << path << ": " << e.what() << '\n';
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
		// must contain LIMA_data or automatedtests
		if (dir.string().find("LIMA_data") == std::string::npos && dir.string().find("automatedtests") == std::string::npos) {
			throw std::runtime_error("LIMA is not allowed to clean this directory");
		}

		if (!fs::exists(dir) || !fs::is_directory(dir)) {
			std::cerr << "Path does not exist or is not a directory: " << dir << std::endl;
			return;
		}

		// Remove all files that do not contain "reference" in their name
		for (auto& p : fs::recursive_directory_iterator(dir)) {
			if (fs::is_regular_file(p) && p.path().filename().string().find(except) == std::string::npos) {
				TryDeleteFile(p.path());
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
	static std::unique_ptr<Environment> basicSetup(const std::string& foldername, std::optional<SimParams> simparams, EnvMode envmode) {
		
		const fs::path work_folder = AutomatedTestsDir() / foldername;
		const GroFile conf{getMostSuitableGroFile(work_folder)};
		const TopologyFile topol {work_folder / "molecule/topol.top"};
		const fs::path simpar = work_folder / "sim_params.txt";

		auto env = std::make_unique<Environment>(work_folder, envmode);

		const SimParams ip = simparams.value_or(SimParams{ simpar });


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

	struct VarianceCoefficientThresholds {
		float max_vc;
		float max_gradient;
	};

	const fs::path VarianceCoefficientTargetsPath() {
		return AutomatedTestsDir() / "vc_targets.csv";
	}

	const fs::path VarianceCoefficientActualsPath() {
		return AutomatedTestsDir() / "vc_results.csv";
	}

	std::map<std::string, VarianceCoefficientThresholds>& ActualVarianceCoefficientResults() {
		static std::map<std::string, VarianceCoefficientThresholds> results;
		return results;
	}

	bool& VarianceCoefficientSuiteStarted() {
		static bool started = false;
		return started;
	}

	void WriteActualVarianceCoefficientResults() {
		std::ofstream output{ VarianceCoefficientActualsPath(), std::ios::trunc };
		if (!output) {
			throw std::runtime_error("Could not write " + VarianceCoefficientActualsPath().string());
		}

		output << "test,max_vc,max_gradient\n" << std::setprecision(9);
		for (const auto& [test_name, result] : ActualVarianceCoefficientResults()) {
			output << test_name << ',' << result.max_vc << ',' << result.max_gradient << '\n';
		}
	}

	void BeginVarianceCoefficientTestSuite() {
		ActualVarianceCoefficientResults().clear();
		VarianceCoefficientSuiteStarted() = true;
		WriteActualVarianceCoefficientResults();
	}

	void RecordActualVarianceCoefficientResult(
		const std::string& test_name,
		const std::vector<float>& VCs,
		const std::vector<float>& energy_gradients
	) {
		if (VCs.empty() || energy_gradients.empty()) {
			throw std::runtime_error("Variance coefficient tests must provide at least one VC and energy gradient");
		}
		if (!VarianceCoefficientSuiteStarted()) {
			BeginVarianceCoefficientTestSuite();
		}

		const float max_vc = *std::max_element(VCs.begin(), VCs.end());
		const float max_gradient = std::abs(*std::max_element(
			energy_gradients.begin(), energy_gradients.end(),
			[](float lhs, float rhs) { return std::abs(lhs) < std::abs(rhs); }
		));

		auto [entry, inserted] = ActualVarianceCoefficientResults().try_emplace(
			test_name, VarianceCoefficientThresholds{ max_vc, max_gradient });
		if (!inserted) {
			entry->second.max_vc = std::max(entry->second.max_vc, max_vc);
			entry->second.max_gradient = std::max(entry->second.max_gradient, max_gradient);
		}
		WriteActualVarianceCoefficientResults();
	}

	const std::map<std::string, VarianceCoefficientThresholds>& VarianceCoefficientTargets() {
		static const auto targets = []() {
			std::ifstream input{ VarianceCoefficientTargetsPath() };
			if (!input) {
				throw std::runtime_error("Could not read " + VarianceCoefficientTargetsPath().string());
			}

			std::string line;
			std::getline(input, line);
			if (line.ends_with('\r')) line.pop_back();
			if (line != "test,max_vc,max_gradient") {
				throw std::runtime_error("Unexpected header in " + VarianceCoefficientTargetsPath().string());
			}

			std::map<std::string, VarianceCoefficientThresholds> parsed;
			while (std::getline(input, line)) {
				if (line.ends_with('\r')) line.pop_back();
				if (line.empty()) continue;

				std::istringstream row{ line };
				std::string test_name, max_vc_text, max_gradient_text;
				if (!std::getline(row, test_name, ',') ||
					!std::getline(row, max_vc_text, ',') ||
					!std::getline(row, max_gradient_text) ||
					test_name.empty()) {
					throw std::runtime_error("Malformed row in " + VarianceCoefficientTargetsPath().string() + ": " + line);
				}

				const VarianceCoefficientThresholds values{ std::stof(max_vc_text), std::stof(max_gradient_text) };
				if (!parsed.emplace(test_name, values).second) {
					throw std::runtime_error("Duplicate variance coefficient target: " + test_name);
				}
			}
			return parsed;
		}();
		return targets;
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
	std::pair<bool, std::string> evaluateTest(
		const std::string& test_name,
		std::vector<float> VCs,
		std::vector<float> energy_gradients)
	{
		RecordActualVarianceCoefficientResult(test_name, VCs, energy_gradients);
		const auto target = VarianceCoefficientTargets().find(test_name);
		if (target == VarianceCoefficientTargets().end()) {
			throw std::runtime_error("No variance coefficient target configured for test: " + test_name);
		}
		const float target_vc = target->second.max_vc;
		const float max_energygradient_abs = target->second.max_gradient;

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
				return { false, std::format("Var. Coeff. of {:.3e} was too far from the target {:.3e}", vc, target_vc) };
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
	static void setConsoleTextColorYellow() { std::cout << "\033[33m"; }
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


		void printStatus(std::string insert="") const {
			if (success) {
				setConsoleTextColorGreen();
			}
			else {
				setConsoleTextColorRed();
			}


			if (error_description.length() > 55) { std::cout << "\n\t"; }
			std::cout << error_description << insert << "\n";


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
				TimeIt timer{};
				std::cout << "Test " << name << " ";
				testresult = std::make_unique<LimaUnittestResult>(test());

				int str_len = 6 + name.length();
				while (str_len++ < 61) { std::cout << " "; }

				//testresult->printStatus(" (" + timer.ElapsedPretty() + ")");
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
		LimaUnittestManager(){ BeginVarianceCoefficientTestSuite(); }
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
		const std::string& test_name,
		std::optional<SimParams> ip = {}
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
		ASSERT(env->getSimPtr()->getStep() == env->getSimPtr()->simParams.n_steps, std::format("Simulation did not finish {}/{}",
			env->getSimPtr()->getStep(), env->getSimPtr()->simParams.n_steps));

		const auto result = evaluateTest(test_name, { varcoff }, {analytics.energy_gradient});

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

	LimaUnittestResult TestIsDeterministic(std::function<LimaUnittestResult()> testFunc, int repetitions, EnvMode envmode) {
		std::vector<std::string> results;
		for (int i = 0; i < repetitions; i++) {
			LimaUnittestResult result = testFunc();
			results.push_back(result.error_description);
		}
		bool allSame = std::all_of(results.begin(), results.end(), [&](const std::string& res) {
			return res == results[0];
			});
		if (!allSame) {
			std::string errorMsg = "Test produced different results in different runs:\n";
			for (size_t i = 0; i < results.size(); i++) {
				errorMsg += std::format("Run {}: {}\n", i + 1, results[i]);
			}
			return LimaUnittestResult{ false, errorMsg, envmode != Headless };
		}
		else {
			return LimaUnittestResult{ true, "Success", envmode != Headless };
		}
	}

	LimaUnittestResult CompareTopologyFiles(const TopologyFile& newTop, const TopologyFile& refTop, EnvMode envmode) {

		auto EqualUnordered = []<std::ranges::input_range R1, std::ranges::input_range R2>(R1&& a, R2&& b) {
			using T = std::ranges::range_value_t<R1>;

			std::vector<T> va(std::ranges::begin(a), std::ranges::end(a));
			std::vector<T> vb(std::ranges::begin(b), std::ranges::end(b));
			if (va.size() != vb.size())
				return false;

			// Bonded interactions are invariant under complete atom-order reversal.
			auto canonicalize = [](T& interaction) {
				auto reversedIds = interaction.ids;
				std::ranges::reverse(reversedIds);
				if (reversedIds < interaction.ids)
					interaction.ids = reversedIds;
			};
			std::ranges::for_each(va, canonicalize);
			std::ranges::for_each(vb, canonicalize);

			auto byIdsAndFunction = [](const T& lhs, const T& rhs) {
				if (lhs.ids != rhs.ids)
					return lhs.ids < rhs.ids;
				return lhs.funct < rhs.funct;
			};
			std::ranges::sort(va, byIdsAndFunction);
			std::ranges::sort(vb, byIdsAndFunction);

			return va == vb;
		};

		ASSERT(std::ranges::equal(newTop.GetAllElements<TopologyFile::AtomsEntry>(), refTop.GetAllElements<TopologyFile::AtomsEntry>()), "Topology AtomsEntry Mismatch");
		ASSERT(EqualUnordered(newTop.GetAllElements<TopologyFile::SingleBond>(), refTop.GetAllElements<TopologyFile::SingleBond>()), "Topology SingleBond Mismatch");
		ASSERT(EqualUnordered(newTop.GetAllElements<TopologyFile::PairBond>(), refTop.GetAllElements<TopologyFile::PairBond>()), "Topology PairBond Mismatch");
		ASSERT(EqualUnordered(newTop.GetAllElements<TopologyFile::AngleBond>(), refTop.GetAllElements<TopologyFile::AngleBond>()), "Topology AngleBond Mismatch");
		ASSERT(EqualUnordered(newTop.GetAllElements<TopologyFile::DihedralBond>(), refTop.GetAllElements<TopologyFile::DihedralBond>()), "Topology DihedralBond Mismatch");
		ASSERT(EqualUnordered(newTop.GetAllElements<TopologyFile::ImproperDihedralBond>(), refTop.GetAllElements<TopologyFile::ImproperDihedralBond>()), "Topology ImproperDihedralBond Mismatch");
		return LimaUnittestResult{ true, "Success", false };
	}

	LimaUnittestResult CompareGroFiles(const GroFile& newGro, const GroFile& refGro, EnvMode envmode, float maxCoordinateError=0.0015) {		
		ASSERT(newGro.box_size == refGro.box_size, "Box size mismatch");
		ASSERT(newGro.atoms.size() == refGro.atoms.size(), "Atom count mismatch");
		for (int i = 0; i < newGro.atoms.size(); i++) {
			const auto& newAtom = newGro.atoms[i];
			const auto& refAtom = refGro.atoms[i];
			ASSERT(newAtom.residue_number == refAtom.residue_number, "Residue number mismatch");
			ASSERT(newAtom.residueName == refAtom.residueName, "Residue name mismatch");
			ASSERT(newAtom.atomName== refAtom.atomName, "Atom name mismatch");
			ASSERT(newAtom.gro_id== refAtom.gro_id, "Atom number mismatch");


			bool errX = std::abs(newAtom.position.x - refAtom.position.x) > maxCoordinateError;
			bool errY = std::abs(newAtom.position.y - refAtom.position.y) > maxCoordinateError;
			bool errZ = std::abs(newAtom.position.z - refAtom.position.z) > maxCoordinateError;
			if (errX || errY || errZ) {
				std::string errorMsg = std::format("Atom {} coordinate mismatch: new ({:.6f}, {:.6f}, {:.6f}) vs ref ({:.6f}, {:.6f}, {:.6f})",
					newAtom.gro_id,
					newAtom.position.x, newAtom.position.y, newAtom.position.z,
					refAtom.position.x, refAtom.position.y, refAtom.position.z);
				return LimaUnittestResult{ false, errorMsg, envmode != Headless };
			}

		}
		return LimaUnittestResult{ true, "Success", false };
	}

} // namespace TestUtils


