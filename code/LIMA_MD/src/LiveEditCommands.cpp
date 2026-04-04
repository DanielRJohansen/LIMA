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





	Command ParseCommand(const std::string& commandString) {
		std::stringstream commandStream(commandString);
		if (commandString.starts_with("insertmolecule ")) {
			return ParseInsertMolecule(commandString);
		}

		

		return Invalid{};
	}



}