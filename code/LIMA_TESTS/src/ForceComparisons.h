#include "TestUtils.h"
#include "Filehandling.h"




namespace ForceComparisons {
	using namespace TestUtils;


	const float Error(const Float3& value, const Float3& ref) {
		return std::abs((value - ref).len() / ref.len());
	}

	std::unique_ptr<Environment> LoadAndRunSim(const fs::path& workdir) {
		GroFile grofile{ workdir / "conf.gro" };
		TopologyFile topfile{ workdir / "topol.top" };
		const fs::path paramsPath = fs::exists(workdir / "sim_params.txt") ? workdir / "sim_params.txt" : workdir.parent_path() / "sim_params.txt";
		SimParams simparams{ paramsPath };
		auto env = std::make_unique<Environment>( workdir, Headless ); 
		env->CreateSimulation(grofile, topfile, simparams);
		env->run();

		return env;
	}

	bool CompareForces(const fs::path& workdir, float errorThreshold) {
		auto env = LoadAndRunSim(workdir);

		const std::vector<Float3> gromacsForces = FileUtils::ReadCsvAsVectorOfFloat3(workdir / "forces.csv");
		const std::vector<Float3> limaForces = env->GetForces(0);

		if (gromacsForces.size() != limaForces.size()) {
			printf("Lima force count %d Gromacs force count %d\n", limaForces.size(), gromacsForces.size());
			return false;
		}

		for (int i = 0; i < limaForces.size(); i++) {
			if (Error(limaForces[i], gromacsForces[i]) > errorThreshold) {
				limaForces[i].print("limaForce");
				gromacsForces[i].print("gromacsForce");
				return false;
			}
		}
		return true;

	}

	bool PoolNoES() {
		const fs::path workdir = simulations_dir / "forcecomparison" / "PoolNoES";
		return CompareForces(workdir, 0.001);
	} 

	bool PoolES() {
		const fs::path workdir = simulations_dir / "forcecomparison" / "PoolES";
		return CompareForces(workdir, 0.001);
	}

	bool Singlebond() {
		const fs::path workdir = simulations_dir / "forcecomparison" / "Singlebond";
		return CompareForces(workdir, 0.001);
	}

	bool T4RmsdAndRmsf() {
		const fs::path workdir = simulations_dir / "forcecomparison" / "T4Lysozyme";
		auto env = LoadAndRunSim(workdir);

		Trajectory traj = env->WriteSimToTrajectory();
		MDFiles::Dump(traj, workdir / "limaTraj.trr");

		const fs::path scriptPath = FileUtils::GetLimaDir() / "dev" / "PyTools" / "CompareRMSD.py";
		const std::string arguments =
			" " + (workdir / "gromacsTraj.trr").string() +
			" " + (workdir / "limaTraj.trr").string() +			
			" " + (workdir / "conf.gro").string();
		const std::string command = "python " + scriptPath.string() + arguments;
		std::system(command.c_str());

		// Call maybe a python script to compare the rmsd and rmsf?
		return true;
	}


	LimaUnittestResult DoAllForceComparisons(EnvMode envmode) {
		ASSERT(PoolNoES(), "PoolNoES failed");
		//ASSERT(PoolES(), "PoolES failed"); // The gromacs part is wrong here
		ASSERT(Singlebond(), "Singlebond failed");

		//ASSERT(T4RmsdAndRmsf(), "T4RmsdAndRmsf failed");

		return LimaUnittestResult{ true, "", envmode==Full };
	}

}