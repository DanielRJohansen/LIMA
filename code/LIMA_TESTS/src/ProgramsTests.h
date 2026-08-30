#include "../LIMA/include/BuildMembrane.h"
#include "TestUtils.h"


namespace ProgramsTests {

	LimaUnittestResult TestBuildMembrane(EnvMode envmode) {
		const fs::path workDir = HeavyTestsDir() / "etc";

		std::vector<std::string> args = {
			"lima", "buildMembrane",
			"-lipids", "DPPC", "60", "DOPC", "40",
			"-centerz", "3.0",
			"-boxsize", "12.0", "10.0","8.0",
			"-emtol", "50",
			"-working_dir", workDir.string().c_str(),
			"-d"
		};

		std::vector<char*> argv;
		for (auto& arg : args)
			argv.push_back(arg.data());

		int ret = buildMembrane(args.size(), argv.data());

		return LimaUnittestResult{ ret == 0, "Success", envmode == Full };
	}


}