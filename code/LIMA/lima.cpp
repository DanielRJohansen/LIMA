#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <format>
#include <iostream>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <vector>

#include "CliDefinitions.h"
#include "CommandlineUtils.h"
#include "Display.h"
#include "Environment.h"
#include "MDFiles.h"
#include "MoleculeUtils.h"
#include "Programs.h"
#include "SetWindowsRegistry.h"
#include "SimulationBuilder.h"
#include "argparser.h"

namespace fs = std::filesystem;

void PrintTiming(const std::chrono::duration<double> enginetime, double totalNsSimulated) {
    const double wall_time_sec = enginetime.count();

    // Calculate performance metrics
    const double ns_per_day = totalNsSimulated / (wall_time_sec / 86400.0);  // 86400 seconds in a day
    const double hr_per_ns = (wall_time_sec / totalNsSimulated) / 3600.0;    // convert to hours per ns

    // Print time and performance info in the GROMACS-like format
    printf("\n");
    printf("               Wall t (s)\n");
    printf("       Time:    %10.3f\n", wall_time_sec);
    printf("                 (ns/day)    (hour/ns)\n");
    printf("Performance:    %10.3f     %10.3f\n", ns_per_day, hr_per_ns);
}

int Cli::RunMdrun(int argc, char** argv) {
    const std::string helpText{ Cli::MdrunHelp };

    ArgParser parser(helpText);

    fs::path workDir;
    fs::path conf = workDir / "conf.gro";
    fs::path topol = workDir / "topol.top";
    fs::path simpar = workDir / "sim_params.txt";
    fs::path conf_out = workDir / "out.gro";
    fs::path trajOut{};

    bool render = false;
	bool uff = false;

    parser.AddOption({ "--conf", "-c", "-conf" }, false, conf);
	parser.AddOption({ "--topology", "-t", "-topology", "-top" }, false, topol);
	parser.AddOption({ "--simparams", "-s", "-simparams" }, false, simpar);
	parser.AddOption({ "--conf-out", "-conf_out", "-co" }, false, conf_out);
	parser.AddOption({ "--trajectory", "-trr", "-traj" }, false, trajOut);
	parser.AddFlag({ "--display", "-d", "-display" }, [&render]() { render = true; });
    parser.AddFlag({ "--uff", "-uff" }, [&uff]() {uff = true; });
    parser.Parse(argc, argv);

    EnvMode envmode = render ? Full : ConsoleOnly;

    auto env = std::make_unique<Environment>(workDir, envmode);

    const SimParams ip(simpar);
    GroFile grofile{ conf };
    TopologyFile topfile{ topol };

    if ((NodeIndex(grofile.box_size.ToInt3()).toFloat3() - grofile.box_size).len() > 0.0001 || grofile.box_size.x != grofile.box_size.y || grofile.box_size.y != grofile.box_size.z) {
        const int newSize = static_cast<int>(std::ceil((std::max)({
            grofile.box_size.x, grofile.box_size.y, grofile.box_size.z })));
        printf("Boxsize was not an integer, or was not cubic. Setting new boxsize to %d", newSize);
        grofile.box_size = Float3(newSize, newSize, newSize);
    }

    env->CreateSimulation(grofile, topfile, ip);

    const std::chrono::duration<double> enginetime = env->run(); // [s]
    printf("Engine time %f\n", enginetime.count());
    env->WriteBoxCoordinatesToFile(grofile);
    grofile.printToFile(conf_out);

    if (!trajOut.empty()) {
        Trajectory traj = env->WriteSimToTrajectory();
		MDFiles::Dump(traj, trajOut);
    }

    if (uff) {
		fs::path path = workDir / "trajectory.uff";
        env->WriteTrajectoryAsUff(path);
    }

    const double total_ns_simulated = static_cast<double>(ip.n_steps) * ip.dt;
	PrintTiming(enginetime, total_ns_simulated);

    return 0;
}

int Cli::RunBuildMembrane(int argc, char** argv) {
    const std::string helpText{ Cli::BuildMembraneHelp };


    EnvMode envmode = ConsoleOnly;
    fs::path workDir{ std::filesystem::current_path() };

    std::vector<std::pair<std::string, double>> lipids; // {name, percentage}
    std::optional<float> membraneCenterZ = std::nullopt;
    Float3 boxsize{};
    float emtol = 100.f;

    ArgParser argparser(helpText);

	argparser.AddOption({ "--lipids", "-lipids" }, true,
        [&lipids](const std::vector<std::string>& args) {
		if (args.size() % 2 != 0) {
			throw CliError("option '--lipids' expects NAME/PERCENT pairs");
		}
		for (std::size_t i = 0; i < args.size(); i += 2) {
            std::string lipidname = args[i];
            double lipidPercentage = 0;
            try {
                lipidPercentage = std::stod(args[i+1]);
            }
            catch (...) {
				throw CliError(std::format("invalid percentage '{}' for lipid '{}'", args[i + 1], args[i]));
            }
			lipids.emplace_back(lipidname, lipidPercentage);
		}
		}
    );

    argparser.AddOption({ "--center-z", "-c", "-centerz" }, false, membraneCenterZ);
	argparser.AddOption({ "--box-size", "-b", "-boxsize" }, true, boxsize, true);
	argparser.AddOption({ "--em-tolerance", "-emtol", "-tolerance" }, false, emtol);
    argparser.AddOption({ "--working-dir", "-working_dir", "-workdir", "-wd" }, false, workDir);
	argparser.AddFlag({ "--display", "-d", "-display" }, [&envmode]() { envmode = Full; });

    argparser.Parse(argc, argv);

	Lipids::Selection lipidselection;
	for (const auto& lipid : lipids) {
		lipidselection.emplace_back(Lipids::Select{ lipid.first, workDir, lipid.second });
	}

    GroFile grofile;
    grofile.box_size = boxsize;
    grofile.title = "Membrane";
    TopologyFile topfile;
    topfile.SetSystem("Membrane");
    SimulationBuilder::CreateMembrane(grofile, topfile, lipidselection, membraneCenterZ.value_or(boxsize.z/2.f));
    auto sim = Programs::EnergyMinimize(grofile, topfile, true, workDir, envmode, true, emtol);

    grofile.printToFile(workDir / "membrane.gro");
    topfile.printToFile(workDir / "membrane.top");

    auto [step, force] = *std::min_element(sim->maxForceBuffer.begin(), sim->maxForceBuffer.end(),
        [](const std::pair<int64_t, float>& a, const std::pair<int64_t, float>& b) {
            return a.second < b.second;
        }
    );
    std::cout << std::format("buildmembrane finished with a min max-force of {:.3f}\n", force);

	return 0;
}

namespace RenderCli {
	inline std::string Lowercase(std::string value) {
		std::ranges::transform(value, value.begin(), [](const unsigned char c) {
			return static_cast<char>(std::tolower(c));
		});
		return value;
	}
}


int Cli::RunRender(int argc, char** argv) {
    const std::string helpText{ Cli::RenderHelp };


    ArgParser parser(helpText);

	fs::path conf = "./conf.gro";
	fs::path topol = "./topol.top";
	std::string water = "tip3p";
    bool whole = false;
    bool hidewater = false;
    std::vector<int> highlightAtomsInput{};

    parser.AddOption({ "--structure", "-f", "--conf", "-c", "-structure", "-conf" }, false, conf);
    parser.AddOption({ "--topology", "-t", "-topology", "-top" }, false, topol);
	parser.AddOption({ "--water-model", "-water-model", "-water" }, false, water);
	parser.AddFlag({ "--whole", "-w", "-whole" }, [&whole]() { whole = true; });
	parser.AddFlag({ "--hide-water", "-hidewater", "-hw" }, [&hidewater]() { hidewater = true; });
	parser.AddOption({ "--highlight", "-highlight", "-hl" }, false, highlightAtomsInput);
    parser.Parse(argc, argv);

	std::optional<Programs::GmxConversionResult> conversion;
	const std::string extension = RenderCli::Lowercase(conf.extension().string());
	if (extension == ".pdb" || extension == ".cif") {
		conversion.emplace(Programs::ToGmx(conf, Programs::ParseWaterModel(water)));
	}
	else if (extension != ".gro") {
		throw std::runtime_error(std::format(
			"lima render expects a .gro, .pdb, or .cif input file, got {}", conf.string()));
	}

	GroFile grofile = conversion ? std::move(conversion->grofile) : GroFile{ conf };

    if (whole) {
		if (conversion) MoleculeUtils::MakeMoleculeWholeAfterPBCFragmentation(grofile, conversion->topology);
		else {
			TopologyFile topfile{ topol };
			MoleculeUtils::MakeMoleculeWholeAfterPBCFragmentation(grofile, topfile);
		}
    }

	// Converted coordinate files carry crystallographic cell dimensions, which do
	// not necessarily bound the displayed biological structure. A whole structure
	// can likewise extend beyond its former periodic cell after unwrapping.
	if (conversion || whole) MoleculeUtils::FitMoleculeInBox(grofile);

	std::set<int> highlightedAtoms(highlightAtomsInput.begin(), highlightAtomsInput.end());

    Display d{};
	std::optional<TopologyFile> topologyFile;
	const TopologyFile* topology = nullptr;
	if (conversion) {
		topology = &conversion->topology;
	}
	else if (fs::exists(topol)) {
		try {
			topologyFile.emplace(topol);
			topology = &topologyFile.value();
		}
		catch (const std::exception&) {
			// An unreadable optional topology should not prevent coordinate rendering.
		}
	}

	std::unique_ptr<Rendering::AtomRenderTask> renderTask;
	if (topology) {
		try {
			// Simulation boxes currently require integer dimensions. Keep the input
			// coordinates unchanged and use a private, padded copy for the adapter.
			GroFile simulationGrofile = grofile;
			simulationGrofile.box_size = Float3{
				std::ceil(simulationGrofile.box_size.x),
				std::ceil(simulationGrofile.box_size.y),
				std::ceil(simulationGrofile.box_size.z)
			};

			Environment environment(simulationGrofile.m_path.parent_path(), EnvMode::Headless);
			environment.CreateSimulation(simulationGrofile, *topology, SimParams{});
			std::unique_ptr<Simulation> simulation = environment.GetSim();
			renderTask = std::make_unique<Rendering::AtomRenderTask>(
				simulation->box->persistentClusters,
				simulation->box->persistentClustersMetadata,
				simulation->box->boxparams,
				SimStatus{}, simulation->box->backboneChains);
		}
		catch (const std::exception&) {
			// Coordinate-only rendering remains valid when the topology or
			// force-field cannot produce a complete simulation.
			renderTask.reset();
		}
	}
	if (!renderTask) {
		renderTask = std::make_unique<Rendering::AtomRenderTask>(grofile, !hidewater);
	}
	renderTask->showSolvents = !hidewater;
	renderTask->highlightedAtoms = highlightedAtoms;
    d.Render(std::move(renderTask), true);

    return 0;
}

int Cli::RunMakeBox(int argc, char** argv) {
    const std::string helpText{ Cli::MakeBoxHelp };


    ArgParser parser(helpText);

    std::string name ="";
    int boxsize{};

    // Add options to parser
    parser.AddOption({ "--box-size", "-b", "-boxsize" }, true, boxsize);
    parser.AddOption({ "--name", "-n", "-name" }, false, name);

    parser.Parse(argc, argv);

	const fs::path confPath = fs::absolute(!name.empty() ? fs::path(name + ".gro") : "conf.gro");
	const fs::path topPath  = fs::absolute(!name.empty() ? fs::path(name + ".top") : "topol.top");

    //printf("Creating files: \n\t%s\n\t%s\n", confPath.string(), topPath.string());


    // Use the paths in your program
    GroFile grofile{};
    grofile.m_path = confPath;
    grofile.box_size = Float3{ static_cast<float>(boxsize) };
    grofile.printToFile();

    TopologyFile topfile{};
    topfile.path = topPath;
    topfile.printToFile();

    return 0;
}

int Cli::RunInsertMolecule(int argc, char** argv) {
    namespace fs = std::filesystem;
    const std::string helpText{ Cli::InsertMoleculeHelp };


    ArgParser parser(helpText);

    fs::path confSrcPath{};
    fs::path topSrcPath{};
    fs::path confTgtPath{"./conf.gro"};
    fs::path topTgtPath{"./topol.top"};

    Float3 position{FLT_MAX, FLT_MAX , FLT_MAX };

    parser.AddOption({ "--conf-source", "-conf_source", "-cs" }, true, confSrcPath);
    parser.AddOption({ "--top-source", "-top_source", "-ts" }, true, topSrcPath);
    parser.AddOption({ "--conf-target", "-conf_target", "-ct" }, false, confTgtPath);
    parser.AddOption({ "--top-target", "-top_target", "-tt" }, false, topTgtPath);
    parser.AddOption({ "--position", "-p", "-position" }, false, position);
    parser.Parse(argc, argv);

    GroFile groSrc{ confSrcPath };
    auto topSrc = std::make_shared<TopologyFile>(topSrcPath);
    GroFile groTgt{ confTgtPath };
    TopologyFile topTgt{ topTgtPath };


    if (position.x == FLT_MAX) {
        position = groTgt.box_size / 2.f;
	}

    SimulationBuilder::InsertSubmoleculeInSimulation(groTgt, topTgt, groSrc, topSrc, position);

    return 0;
}

int Cli::RunInsertMolecules(int argc, char** argv) {
    namespace fs = std::filesystem;
    const std::string helpText{ Cli::InsertMoleculesHelp };


    ArgParser parser(helpText);

    fs::path confSrcPath{};
    fs::path topSrcPath{};
    fs::path confTgtPath{"./conf.gro"};
    fs::path topTgtPath{"./topol.top"};

    bool rotateRandomly = false;
    int nInsertions{};
    bool display = false;

    parser.AddOption({ "--conf-source", "-conf_source", "-cs" }, true, confSrcPath);
    parser.AddOption({ "--top-source", "-top_source", "-ts" }, true, topSrcPath);
    parser.AddOption({ "--conf-target", "-conf_target", "-ct" }, false, confTgtPath);
    parser.AddOption({ "--top-target", "-top_target", "-tt" }, false, topTgtPath);
    parser.AddOption({ "--num-insertions", "-n", "-num_insertions" }, true, nInsertions);
    parser.AddFlag({ "--rotate-randomly", "-rotate_randomly", "-rr" }, [&rotateRandomly]() {rotateRandomly = true; });
    parser.AddFlag({ "--display", "-d", "-display" }, [&display]() {display = true; });
    parser.Parse(argc, argv);

    if (!fs::exists(confSrcPath)) throw std::runtime_error(std::format("source coordinates not found: {}", confSrcPath.string()));
    if (!fs::exists(topSrcPath)) throw std::runtime_error(std::format("source topology not found: {}", topSrcPath.string()));
    if (!fs::exists(confTgtPath)) throw std::runtime_error(std::format("target coordinates not found: {}", confTgtPath.string()));
	if (!fs::exists(topTgtPath)) throw std::runtime_error(std::format("target topology not found: {}", topTgtPath.string()));

	confTgtPath = fs::absolute(confTgtPath);
	topTgtPath = fs::absolute(topTgtPath);

    GroFile groSrc{ confSrcPath };
    auto topSrc = std::make_shared<TopologyFile>(topSrcPath);
    GroFile groTgt{ confTgtPath };
    TopologyFile topTgt{};
    topTgt.SetSystem(topSrc->GetSystem().title + " " + std::to_string(nInsertions));


    Programs::EnergyMinimize(groSrc, *topSrc, true, fs::current_path(), display ? Full : Headless, false);

    SimulationBuilder::InsertSubmoleculesInSimulation(groTgt, topTgt, groSrc, topSrc, nInsertions, rotateRandomly);
    Programs::StaticbodyEnergyMinimize(groTgt, topTgt, display);
    groTgt.printToFile(confTgtPath);
    topTgt.printToFile(topTgtPath);

    return 0;
}

int Cli::RunEditConf(int argc, char** argv) {
    const std::string helpText{ Cli::EditConfHelp };


    ArgParser parser(helpText);

    fs::path confInputPath{};
    fs::path topInputPath{};
    fs::path confOutputPath{};

    bool makeWhole = false;
    Float3 setCenter{FLT_MAX, FLT_MAX , FLT_MAX };
    Float3 rotate{FLT_MAX, FLT_MAX , FLT_MAX };

    parser.AddOption({ "--conf", "-c", "-conf" }, true, confInputPath);
    parser.AddOption({ "--topology", "-t", "-top" }, true, topInputPath);
    parser.AddOption({ "--conf-out", "-conf_out", "-co" }, false, confOutputPath);
    parser.AddFlag({ "--whole", "-w", "-whole" }, [&makeWhole]() {makeWhole = true; });
    parser.AddOption({ "--set-center", "-set_center", "-sc" }, false, setCenter);
    parser.AddOption({ "--rotate", "-r", "-rotate" }, false, rotate);
    parser.Parse(argc, argv);

    if (confOutputPath.empty())
        confOutputPath = confInputPath;


    GroFile grofile{ confInputPath };
    TopologyFile topfile{ topInputPath };

    bool setCenterChosen = setCenter.x != FLT_MAX;
    bool rotateChosen = rotate.x != FLT_MAX;

    if (makeWhole || setCenterChosen || rotateChosen) {
        MoleculeUtils::MakeMoleculeWholeAfterPBCFragmentation(grofile, topfile.GetMoleculeType());
    }

    if (setCenterChosen) {
        MoleculeUtils::CenterMolecule(grofile, topfile.GetMoleculeType(), setCenter);
    }

    if (rotateChosen) {
        MoleculeUtils::RotateMolecule(grofile, rotate);
	}

    grofile.printToFile(confOutputPath);

    return 0;
}

int Cli::RunEnergyMinimization(int argc, char** argv) {
    namespace fs = std::filesystem;
    const std::string helpText{ Cli::EnergyMinimizationHelp };


    ArgParser parser(helpText);

    fs::path confPath{"conf.gro"};
    fs::path topPath{"topol.top"};
    fs::path confPathOut{"conf.gro"};

    bool render = false;
    float emtol = 100.f;

    parser.AddOption({ "--conf", "-c", "-conf" }, false, confPath);
    parser.AddOption({ "--topology", "-t", "-top" }, false, topPath);
    parser.AddOption({ "--conf-out", "-conf_output", "-co" }, false, confPathOut);
    parser.AddOption({ "--em-tolerance", "-emtol" }, false, emtol);
    parser.AddFlag({ "--display", "-d", "-display" }, [&render](){render=true;});
    parser.Parse(argc, argv);

    GroFile grofile{ fs::canonical(confPath) };
    TopologyFile topfile{ fs::canonical(topPath) };

    Programs::EnergyMinimize(grofile, topfile, true, fs::current_path(), render ? Full : ConsoleOnly, false, emtol);

    grofile.printToFile(fs::absolute(confPathOut));

    return 0;
}

int Cli::RunToGmx(int argc, char** argv) {
    const std::string helpText{ Cli::ToGmxHelp };

    ArgParser parser(helpText);
    fs::path structurePath;
    std::string name;
    std::string water = "tip3p";
    parser.AddOption({ "--structure", "-f", "-structure" }, true, structurePath);
    parser.AddOption({ "--name", "-n", "-name" }, false, name);
    parser.AddOption({ "--water-model", "-water-model", "-water" }, false, water);
    parser.Parse(argc, argv);

    if (!name.empty() && fs::path{ name }.has_parent_path()) {
        throw std::runtime_error("togmx output name must be a basename");
    }

    auto conversion = Programs::ToGmx(structurePath, Programs::ParseWaterModel(water));
    const fs::path directory = structurePath.parent_path().empty() ? fs::current_path() : structurePath.parent_path();
    const fs::path groPath = directory / (name.empty() ? "conf.gro" : name + ".gro");
    const fs::path topPath = directory / (name.empty() ? "topol.top" : name + ".top");

    const auto& molecules = conversion.topology.GetSystem().molecules;
    if (molecules.size() != conversion.positionRestraints.size()) {
        throw std::runtime_error("Converted molecule and position-restraint counts differ");
    }
    for (std::size_t i = 0; i < molecules.size(); ++i) {
        const std::string suffix = i == 0 ? std::string{} : "_" + molecules[i].name;
        const fs::path filename = name.empty()
            ? fs::path{ "posre" + suffix + ".itp" }
            : fs::path{ name + "_posre" + suffix + ".itp" };
        conversion.topology.moleculetypes.at(molecules[i].name)->positionRestraintsInclude = filename;
        conversion.positionRestraints[i].printToFile(directory / filename);
    }

    conversion.grofile.printToFile(groPath);
    conversion.topology.printToFile(topPath);
    return 0;
}

void SelfTest() {
	const std::filesystem::path workDir = std::filesystem::current_path() / "selftest";

	const fs::path slipidsPath = FileUtils::GetLimaDir() / "resources/Slipids";
	std::vector<std::string> targets;
	for (const auto& entry : fs::directory_iterator(slipidsPath)) {
		if (entry.path().extension() == ".gro") {
			std::string base_name = entry.path().stem().string();
			if (fs::exists(slipidsPath / (base_name + ".itp"))) {
				targets.push_back(base_name);
			}
		}
	}

	Lipids::Selection lipidselection;
	for (const auto& lipidname : targets) {
		lipidselection.emplace_back(Lipids::Select{ lipidname, workDir, 100. / static_cast<double>(targets.size()) });
	}

	GroFile gro;
	gro.box_size = Float3{ 10.f };
	gro.title = "Membrane";
	TopologyFile top;
	top.SetSystem("Membrane");
	SimulationBuilder::CreateMembrane(gro, top, lipidselection, 5.f);
	Programs::EnergyMinimize(gro, top, false, workDir, Full, true, 5000.f);

	printf("Selftest successful"); // Otherwise we'd have thrown by now
}

int Cli::RunMakeSimParams(int argc, char** argv) {
    ArgParser parser{ std::string{ MakeSimParamsHelp } };
    parser.Parse(argc, argv);
    SimParams params{};
    params.DumpToFile();
    return 0;
}

int Cli::RunSelfTest(int argc, char** argv) {
    ArgParser parser{ std::string{ SelfTestHelp } };
    parser.Parse(argc, argv);
    SelfTest();
    return 0;
}

namespace {

const Cli::CommandDefinition* FindCommand(const std::string_view name) {
    const auto command = std::ranges::find(Cli::Commands, name, &Cli::CommandDefinition::name);
    return command == Cli::Commands.end() ? nullptr : &*command;
}

void PrintGeneralHelp() {
    std::cout << "Usage: lima COMMAND [OPTION]...\n"
                 "       lima [--help] [--version]\n\n"
                 "LIMA is a suite of molecular-dynamics and membrane-simulation tools.\n\n"
                 "Commands:\n";
    for (const auto& command : Cli::Commands)
        std::cout << std::format("  {:<20} {}\n", command.name, command.summary);
    std::cout << "\nRun 'lima help COMMAND' or 'lima COMMAND --help' for command help.\n";
}

int Dispatch(int argc, char** argv) {
    if (argc == 1) {
        PrintGeneralHelp();
        return 0;
    }

    const std::string_view argument{ argv[1] };
    if (argument == "--help" || argument == "-h" || argument == "-help") {
        PrintGeneralHelp();
        return 0;
    }
    if (argument == "--version") {
        std::cout << "lima (development build)\n";
        return 0;
    }
    if (argument == "help") {
        if (argc == 2) {
            PrintGeneralHelp();
            return 0;
        }
        if (argc != 3) throw CliError("usage: lima help COMMAND");
        const auto* command = FindCommand(argv[2]);
        if (!command) throw CliError(std::format("unrecognized command '{}'", argv[2]));
        std::cout << command->helpText;
        return 0;
    }
    if (argument == "setregistry") {
        if (argc != 2) throw CliError("'setregistry' does not accept options");
        RegisterGrofileAssociation();
        return 0;
    }

    const auto* command = FindCommand(argument);
    if (!command) throw CliError(std::format("unrecognized command '{}'", argument));
    return command->handler(argc, argv);
}

} // namespace

int main(int argc, char** argv) {
    try {
        return Dispatch(argc, argv);
    }
    catch (const HelpRequested& help) {
        std::cout << help.helpText;
        return 0;
    }
    catch (const CliError& error) {
        std::cerr << "lima: " << error.what() << "\nTry 'lima --help' for more information.\n";
        return 2;
    }
    catch (const std::exception& error) {
        std::cerr << "lima: " << error.what() << '\n';
        return 1;
    }
    catch (...) {
        std::cerr << "lima: unknown error\n";
        return 1;
    }
}
