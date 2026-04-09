#include "LiveEditCommands.h"
#include "argparser.h"


namespace LiveEdit {





	InsertMolecule ParseInsertMolecule(const std::string& _commandStr) {
		//std::string word

		std::vector<std::string> args;
		args.push_back("lima"); // argsparser exptects first arg to be the program, so we emulate that behaviour here..
		std::stringstream commandStream(_commandStr);
		std::string word;
		while (commandStream >> word) {
			args.push_back(word);
		}

		ArgParser parser{""};
		fs::path conf = "./conf.gro";
		fs::path topol = "./topol.top";
		Float3 _pos(FLT_MAX, FLT_MAX, FLT_MAX);

		parser.AddOption({ "-conf", "-c" }, false, conf);
		parser.AddOption({ "-topology", "-top", "-t" }, false, topol);
		parser.AddOption({ "-position", "-p" }, false, _pos, true);
		parser.Parse(args);

		std::optional<Float3> pos = _pos.x == FLT_MAX ? std::nullopt : std::optional<Float3>{ _pos };		

		return InsertMolecule{ conf, topol, pos };
	}


	BuildMembrane ParseCreateMembrane(const std::string& _commandStr) {
		std::vector<std::string> args;
		args.push_back("lima"); // argsparser exptects first arg to be the program, so we emulate that behaviour here..
		std::stringstream commandStream(_commandStr);
		std::string word;
		while (commandStream >> word) {
			args.push_back(word);
		}
		ArgParser parser{ "" };
		std::vector<std::tuple<std::string, double>> lipids; // {name, percentage}
		std::optional<float> membraneCenterZ = std::nullopt;
		parser.AddOption({ "-lipids" }, true,
			[&lipids](const std::vector<std::string>& args) {
				if (args.size() % 2 != 0) {
					throw std::runtime_error("Invalid -lipids argument. It must have a multiple of two values.");
				}
				for (size_t i = 0; i < args.size(); i += 2) {
					std::string lipidname = args[i];
					double lipidPercentage = 0;
					try {
						lipidPercentage = std::stod(args[i + 1]);
					}
					catch (...) {
						throw std::runtime_error(std::format("Invalid lipid percentage: '{}'. Expected a floating-point value.", args[i]));
						exit(1);
					}
					lipids.emplace_back(lipidname, lipidPercentage);
				}
			}
		);
		parser.AddOption({ "-centerz", "-c" }, false, membraneCenterZ);
		parser.Parse(args);
		return BuildMembrane{ lipids, membraneCenterZ };
	}


	Command ParseCommand(const std::string& commandString) {

		// TODO: Handle capitization for all params, except filenames, since linux cares there.
		std::stringstream commandStream(commandString);
		if (commandString.starts_with("insertmolecule ")) {
			return ParseInsertMolecule(commandString);
		}
		if (commandString.starts_with("buildmembrane ")) {
			return ParseCreateMembrane(commandString);
		}

		return Invalid{};
	}



}