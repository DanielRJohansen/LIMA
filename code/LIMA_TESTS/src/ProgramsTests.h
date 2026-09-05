#include "../LIMA/include/BuildMembrane.h"
#include "Programs.h"
#include "TestUtils.h"

#include <cctype>
#include <set>

namespace ProgramsTests {
	namespace {
		std::set<int> ReadPositionRestraintAtoms(const fs::path& path) {
			std::ifstream input(path);
			if (!input) throw std::runtime_error(std::format("Could not open position restraints {}", path.string()));
			std::set<int> result;
			for (std::string line; std::getline(input, line);) {
				if (const auto comment = line.find(';'); comment != std::string::npos) line.erase(comment);
				const auto fields = StringUtils::SplitWords(line);
				if (fields.size() == 5 && std::isdigit(static_cast<unsigned char>(fields[0][0]))) result.insert(std::stoi(fields[0]));
			}
			return result;
		}
	}

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

	LimaUnittestResult TestPdb2Gmx_pdbfile(EnvMode envmode) {
		const fs::path directory = AutomatedTestsDir() / "pdb2gmx";
		const fs::path generatedGro = directory / "generated.gro";
		const fs::path generatedTop = directory / "generated.top";
		const fs::path generatedPosre = directory / "generated_posre.itp";
		TryDeleteFile(generatedGro);
		TryDeleteFile(generatedTop);
		TryDeleteFile(generatedPosre);

		Programs::pdb2gmx(directory / "6lzm.pdb", "generated");
		LimaUnittestResult gResult = CompareGroFiles(GroFile{ generatedGro }, GroFile{ directory / "conf_ref.gro" }, envmode);
		if (!gResult.success) 
			return LimaUnittestResult{ false, gResult.error_description, envmode == Full };
		
		LimaUnittestResult tResult = CompareTopologyFiles(generatedTop, directory / "topol_ref.top", envmode);	
		if (!tResult.success) 
			return LimaUnittestResult{ false, tResult.error_description, envmode == Full };

		ASSERT(ReadPositionRestraintAtoms(generatedPosre) == ReadPositionRestraintAtoms(directory / "posre_ref.itp"),
			"Position-restraint atom sets differ");
		return LimaUnittestResult{ true, "Success", envmode == Full };
	}
}
