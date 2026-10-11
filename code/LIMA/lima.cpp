#include <algorithm>
#include "Format.h"
#include <cctype>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <vector>

#include "CliDefinitions.h"
#include "Daemon.h"
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

void PrintTiming(std::ostream& out, const std::chrono::duration<double> enginetime, double totalNsSimulated) {
    const double wall_time_sec = enginetime.count();

    // Calculate performance metrics
    const double ns_per_day = totalNsSimulated / (wall_time_sec / 86400.0);  // 86400 seconds in a day
    const double hr_per_ns = (wall_time_sec / totalNsSimulated) / 3600.0;    // convert to hours per ns

    // Print time and performance info in the GROMACS-like format
    out << "\n";
    out << "               Wall t (s)\n";
    out << Lima::Format("       Time:    {:10.3f}\n", wall_time_sec);
    out << "                 (ns/day)    (hour/ns)\n";
    out << Lima::Format("Performance:    {:10.3f}     {:10.3f}\n", ns_per_day, hr_per_ns);
}

int Cli::RunMdrun(const Context& ctx, int argc, char** argv) {
    const std::string helpText{ Cli::MdrunHelp };

    ArgParser parser(helpText);

    const fs::path& workDir = ctx.cwd;
    fs::path conf = "conf.gro";
    fs::path topol = "topol.top";
    fs::path simpar = "sim_params.txt";
    fs::path conf_out = "out.gro";
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
    conf = ctx.Resolve(conf);
    topol = ctx.Resolve(topol);
    simpar = ctx.Resolve(simpar);
    conf_out = ctx.Resolve(conf_out);
    trajOut = ctx.Resolve(trajOut);

    EnvMode envmode = render ? Full : ConsoleOnly;

    Environment& env = Environment::Get();

    const SimParams ip(simpar);

    // The files are parsed by the Environment, so this command joins a batch as early as possible
    SimulationJob job;
    job.workDir = workDir;
    job.grofile.emplace(conf);
    job.topfile.emplace(topol);
    job.simParams = ip;
    job.mode = envmode;
    std::optional<int> resizedBox;
    job.preprocess = [&resizedBox](GroFile& grofile, TopologyFile&, SimParams&) {
        if ((NodeIndex(grofile.box_size.ToInt3()).toFloat3() - grofile.box_size).len() > 0.0001 || grofile.box_size.x != grofile.box_size.y || grofile.box_size.y != grofile.box_size.z) {
            resizedBox = static_cast<int>(std::ceil((std::max)({
                grofile.box_size.x, grofile.box_size.y, grofile.box_size.z })));
            grofile.box_size = Float3(*resizedBox, *resizedBox, *resizedBox);
        }
    };
    auto result = env.Submit(std::move(job)).Get();
    if (resizedBox)
        ctx.out << Lima::Format("Boxsize was not an integer, or was not cubic. Setting new boxsize to {}\n", *resizedBox);

    const std::chrono::duration<double> enginetime = result.engineTime; // [s]
    ctx.out << Lima::Format("Engine time {:f}\n", enginetime.count());
    if (result.execution.batchSize > 1)
        ctx.out << Lima::Format("Ran in a batch of {} simulations\n", result.execution.batchSize);
    result.FinalSystem().coordinates.printToFile(conf_out);

    if (!trajOut.empty()) {
        Trajectory traj = result.MakeTrajectory();
		MDFiles::Dump(traj, trajOut);
    }

    if (uff) {
		fs::path path = workDir / "trajectory.uff";
        result.WriteTrajectoryAsUff(path);
    }

    const double total_ns_simulated = static_cast<double>(ip.n_steps) * ip.dt;
	PrintTiming(ctx.out, enginetime, total_ns_simulated);

    return 0;
}

int Cli::RunBuildMembrane(const Context& ctx, int argc, char** argv) {
    const std::string helpText{ Cli::BuildMembraneHelp };


    EnvMode envmode = ConsoleOnly;
    fs::path workDir{ ctx.cwd };

    std::vector<std::pair<std::string, double>> lipids; // {name, percentage}
    std::optional<float> membraneCenterZ = std::nullopt;
    std::optional<MembraneGeometry::Sphere> sphere;
    std::optional<MembraneGeometry::Ellipsoid> ellipsoid;
    Float3 boxsize{};
    float emtol = 100.f;
    int randomSeed = 0;

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
				throw CliError(Lima::Format("invalid percentage '{}' for lipid '{}'", args[i + 1], args[i]));
            }
			lipids.emplace_back(lipidname, lipidPercentage);
		}
		}
    );

    argparser.AddOption({ "--plane", "--center-z", "-c", "-centerz", "-plane" }, false, membraneCenterZ);
	argparser.AddOption({ "--sphere", "-sphere" }, false,
		[&sphere](const std::vector<std::string>& args) {
			try {
				sphere = MembraneGeometry::Sphere{
					Float3{ std::stof(args[0]), std::stof(args[1]), std::stof(args[2]) },
					std::stof(args[3]) };
			}
			catch (...) { throw CliError("option '--sphere' contains an invalid number"); }
		}, 4, 4);
	argparser.AddOption({ "--ellipsoid", "-ellipsoid" }, false,
		[&ellipsoid](const std::vector<std::string>& args) {
			try {
				ellipsoid = MembraneGeometry::Ellipsoid{
					Float3{ std::stof(args[0]), std::stof(args[1]), std::stof(args[2]) },
					Float3{ std::stof(args[3]), std::stof(args[4]), std::stof(args[5]) } };
			}
			catch (...) { throw CliError("option '--ellipsoid' contains an invalid number"); }
		}, 6, 6);
	argparser.AddOption({ "--box-size", "-b", "-boxsize" }, true, boxsize, true);
	argparser.AddOption({ "--em-tolerance", "-emtol", "-tolerance" }, false, emtol);
	argparser.AddOption({ "--seed", "-seed" }, false, randomSeed);
    argparser.AddOption({ "--working-dir", "-working_dir", "-workdir", "-wd" }, false, workDir);
	argparser.AddFlag({ "--display", "-d", "-display" }, [&envmode]() { envmode = Full; });

    argparser.Parse(argc, argv);
    workDir = ctx.Resolve(workDir);

	const int geometryCount = static_cast<int>(membraneCenterZ.has_value())
		+ static_cast<int>(sphere.has_value()) + static_cast<int>(ellipsoid.has_value());
	if (geometryCount > 1)
		throw CliError("specify only one membrane geometry: --plane, --sphere, or --ellipsoid");
	if (sphere && sphere->radius <= 0.f)
		throw CliError("option '--sphere' requires a positive radius");
	if (ellipsoid && (ellipsoid->radii.x <= 0.f || ellipsoid->radii.y <= 0.f || ellipsoid->radii.z <= 0.f))
		throw CliError("option '--ellipsoid' requires positive radii");

	Lipids::Composition composition{ {}, workDir };
	for (const auto& [name, percentage] : lipids)
		composition.lipids.push_back({ name, percentage });

    MembraneGeometry::Figure geometry = MembraneGeometry::Plane{ membraneCenterZ.value_or(boxsize.z / 2.f) };
    if (sphere) geometry = *sphere;
    else if (ellipsoid) geometry = *ellipsoid;
    // The lipids are loaded and the membrane built by the Environment, so this command joins a batch as early as possible
    auto job = Programs::MakeMembraneJob(workDir, std::move(composition), boxsize, geometry, randomSeed,
        SimParams::BasicEMSimParams(emtol), envmode);
    job.grofile->Get().title = "Membrane";
    job.topfile->Get().SetSystem("Membrane");
    auto emResult = Environment::Get().Submit(std::move(job)).Get();

    const MolecularSystem membrane = emResult.FinalSystem();
    membrane.coordinates.printToFile(workDir / "membrane.gro");
    membrane.topology.printToFile(workDir / "membrane.top");

    auto [step, force] = *std::min_element(emResult.simulation->maxForceBuffer.begin(), emResult.simulation->maxForceBuffer.end(),
        [](const std::pair<int64_t, float>& a, const std::pair<int64_t, float>& b) {
            return a.second < b.second;
        }
    );
    ctx.out << Lima::Format("buildmembrane finished with a min max-force of {:.3f}\n", force);

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


int Cli::RunRender(const Context& ctx, int argc, char** argv) {
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
    conf = ctx.Resolve(conf);
    topol = ctx.Resolve(topol);

	std::optional<Programs::GmxConversionResult> conversion;
	const std::string extension = RenderCli::Lowercase(conf.extension().string());
	if (extension == ".pdb" || extension == ".cif") {
		conversion.emplace(Programs::ToGmx(conf, Programs::ParseWaterModel(water)));
	}
	else if (extension != ".gro") {
		throw std::runtime_error(Lima::Format(
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

			Environment& environment = Environment::Get();
			SimulationJob job;
			job.workDir = simulationGrofile.m_path.parent_path();
			job.grofile = std::move(simulationGrofile);
			job.topfile.emplace(*topology);
			job.simParams = SimParams{};
			job.run = false;
			std::unique_ptr<Simulation> simulation = environment.Submit(std::move(job)).Get().simulation;
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
    d.Submit(1, std::move(renderTask), true);

    return 0;
}

int Cli::RunMakeBox(const Context& ctx, int argc, char** argv) {
    const std::string helpText{ Cli::MakeBoxHelp };


    ArgParser parser(helpText);

    std::string name ="";
    int boxsize{};

    // Add options to parser
    parser.AddOption({ "--box-size", "-b", "-boxsize" }, true, boxsize);
    parser.AddOption({ "--name", "-n", "-name" }, false, name);

    parser.Parse(argc, argv);

	const fs::path confPath = ctx.Resolve(!name.empty() ? fs::path(name + ".gro") : "conf.gro");
	const fs::path topPath  = ctx.Resolve(!name.empty() ? fs::path(name + ".top") : "topol.top");

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

int Cli::RunSolvate(const Context& ctx, int argc, char** argv) {
    ArgParser parser{ std::string{ SolvateHelp } };

    fs::path confPath = "conf.gro";
    fs::path topologyPath = "topol.top";
    fs::path confOutputPath;
    fs::path topologyOutputPath;
    int density = SimulationBuilder::defaultSolventsPerNm3;

    parser.AddOption({ "--conf", "-c", "-conf", "-confInput", "-ci" }, false, confPath);
    parser.AddOption({ "--topology", "-t", "-top", "-topInput", "-ti" }, false, topologyPath);
    parser.AddOption({ "--conf-out", "-co", "-confOutput" }, false, confOutputPath);
    parser.AddOption({ "--topology-out", "-to", "-topOutput" }, false, topologyOutputPath);
    parser.AddOption({ "--density", "-p", "-pressure" }, false, density);
    parser.Parse(argc, argv);
    confPath = ctx.Resolve(confPath);
    topologyPath = ctx.Resolve(topologyPath);
    confOutputPath = ctx.Resolve(confOutputPath);
    topologyOutputPath = ctx.Resolve(topologyOutputPath);

    if (!fs::exists(confPath))
        throw CliError(Lima::Format("input coordinates not found: {}", confPath.string()));
    if (!fs::exists(topologyPath))
        throw CliError(Lima::Format("input topology not found: {}", topologyPath.string()));
    if (density <= 0)
        throw CliError("option '--density' must be positive");

    if (confOutputPath.empty())
        confOutputPath = confPath.parent_path() / (confPath.stem().string() + "_solvated.gro");
    if (topologyOutputPath.empty())
        topologyOutputPath = topologyPath.parent_path() / (topologyPath.stem().string() + "_solvated.top");

    GroFile grofile{ confPath };
    TopologyFile topfile{ topologyPath };
    if (!topfile.HasSystem()) topfile.SetSystem("Solvated system");
    SimulationBuilder::SolvateGrofile(grofile, topfile, density);
    grofile.printToFile(confOutputPath);
    topfile.printToFile(topologyOutputPath);
    return 0;
}

int Cli::RunInsertMolecule(const Context& ctx, int argc, char** argv) {
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
    confSrcPath = ctx.Resolve(confSrcPath);
    topSrcPath = ctx.Resolve(topSrcPath);
    confTgtPath = ctx.Resolve(confTgtPath);
    topTgtPath = ctx.Resolve(topTgtPath);

    GroFile groSrc{ confSrcPath };
    auto topSrc = std::make_shared<TopologyFile>(topSrcPath);
    GroFile groTgt{ confTgtPath };
    TopologyFile topTgt{ topTgtPath };
    // Eg. an empty box from 'lima makebox' has no system yet
    if (!topTgt.HasSystem())
        topTgt.SetSystem(topSrc->HasSystem() ? topSrc->GetSystem().title : "System");

    if (position.x == FLT_MAX) {
        position = groTgt.box_size / 2.f;
	}

    SimulationBuilder::InsertSubmoleculeInSimulation(groTgt, topTgt, groSrc, topSrc, position);
    groTgt.printToFile(confTgtPath);
    topTgt.printToFile(topTgtPath);

    return 0;
}

int Cli::RunInsertMolecules(const Context& ctx, int argc, char** argv) {
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
    confSrcPath = ctx.Resolve(confSrcPath);
    topSrcPath = ctx.Resolve(topSrcPath);
    confTgtPath = ctx.Resolve(confTgtPath);
    topTgtPath = ctx.Resolve(topTgtPath);

    if (!fs::exists(confSrcPath)) throw std::runtime_error(Lima::Format("source coordinates not found: {}", confSrcPath.string()));
    if (!fs::exists(topSrcPath)) throw std::runtime_error(Lima::Format("source topology not found: {}", topSrcPath.string()));
    if (!fs::exists(confTgtPath)) throw std::runtime_error(Lima::Format("target coordinates not found: {}", confTgtPath.string()));
	if (!fs::exists(topTgtPath)) throw std::runtime_error(Lima::Format("target topology not found: {}", topTgtPath.string()));

    GroFile groSrc{ confSrcPath };
    auto topSrc = std::make_shared<TopologyFile>(topSrcPath);
    GroFile groTgt{ confTgtPath };
    // Keep the molecules already in the target, so the topology still matches the coordinates
    TopologyFile topTgt{ topTgtPath };
    if (!topTgt.HasSystem())
        topTgt.SetSystem(topSrc->GetSystem().title + " " + std::to_string(nInsertions));


    // Always headless: a Full-mode job leaves the Environment's Display alive, and StaticbodyEnergyMinimize opens
    // its own. Two Displays drive GLFW from two threads at once, which crashes
    auto emResult = Environment::Get().Submit(SimulationJob{
        ctx.cwd, groSrc, *topSrc, SimParams::BasicEMSimParams(), Headless }).Get();
    emResult.WriteCoordinatesTo(groSrc);

    SimulationBuilder::InsertSubmoleculesInSimulation(groTgt, topTgt, groSrc, topSrc, nInsertions, rotateRandomly);
    Programs::StaticbodyEnergyMinimize(groTgt, topTgt, display);
    groTgt.printToFile(confTgtPath);
    topTgt.printToFile(topTgtPath);

    return 0;
}

int Cli::RunEditConf(const Context& ctx, int argc, char** argv) {
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
    confInputPath = ctx.Resolve(confInputPath);
    topInputPath = ctx.Resolve(topInputPath);
    confOutputPath = ctx.Resolve(confOutputPath);

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

int Cli::RunEnergyMinimization(const Context& ctx, int argc, char** argv) {
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

    auto emResult = Environment::Get().Submit(SimulationJob{ ctx.cwd, fs::canonical(ctx.Resolve(confPath)),
        fs::canonical(ctx.Resolve(topPath)), SimParams::BasicEMSimParams(emtol), render ? Full : ConsoleOnly }).Get();
    emResult.FinalSystem().coordinates.printToFile(ctx.Resolve(confPathOut));

    return 0;
}

int Cli::RunToGmx(const Context& ctx, int argc, char** argv) {
    const std::string helpText{ Cli::ToGmxHelp };

    ArgParser parser(helpText);
    fs::path structurePath;
    std::string name;
    std::string water = "tip3p";
    parser.AddOption({ "--structure", "-f", "-structure" }, true, structurePath);
    parser.AddOption({ "--name", "-n", "-name" }, false, name);
    parser.AddOption({ "--water-model", "-water-model", "-water" }, false, water);
    parser.Parse(argc, argv);
    structurePath = ctx.Resolve(structurePath);

    if (!name.empty() && fs::path{ name }.has_parent_path()) {
        throw std::runtime_error("togmx output name must be a basename");
    }

    auto conversion = Programs::ToGmx(structurePath, Programs::ParseWaterModel(water));
    const fs::path directory = structurePath.parent_path();
    const fs::path groPath = directory / (name.empty() ? "conf.gro" : name + ".gro");
    const fs::path topPath = directory / (name.empty() ? "topol.top" : name + ".top");

    const auto molecules = conversion.topology.GetSystem().Instances() | std::ranges::to<std::vector<TopologyFile::MoleculeEntry>>();
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

void SelfTest(const Cli::Context& ctx) {
	const std::filesystem::path workDir = ctx.cwd / "selftest";

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

	Lipids::Composition composition{ {}, workDir };
	for (const auto& lipidname : targets)
		composition.lipids.push_back({ lipidname, 100. / static_cast<double>(targets.size()) });

	auto job = Programs::MakeMembraneJob(workDir, std::move(composition), Float3{ 10.f }, MembraneGeometry::Plane{ 5.f },
		0, SimParams::BasicEMSimParams(5000.f), Full);
	job.grofile->Get().title = "Membrane";
	job.topfile->Get().SetSystem("Membrane");
	Environment::Get().Submit(std::move(job)).Get();

	ctx.out << "Selftest successful\n"; // Otherwise we'd have thrown by now
}

int Cli::RunMakeSimParams(const Context& ctx, int argc, char** argv) {
    ArgParser parser{ std::string{ MakeSimParamsHelp } };
    parser.Parse(argc, argv);
    SimParams params{};
    params.DumpToFile(ctx.cwd / "sim_params.txt");
    return 0;
}

int Cli::RunSelfTest(const Context& ctx, int argc, char** argv) {
    ArgParser parser{ std::string{ SelfTestHelp } };
    parser.Parse(argc, argv);
    SelfTest(ctx);
    return 0;
}

int Cli::RunWait(const Context& ctx, int argc, char** argv) {
    ArgParser parser{ std::string{ WaitHelp } };
    parser.Parse(argc, argv);
    return Daemon::RunWait(ctx, std::vector<std::string>(argv, argv + argc));
}

int Cli::RunDaemon(const Context& ctx, int argc, char** argv) {
    const std::vector<std::string> args(argv, argv + argc);
    if (args.size() == 3 && (args[2] == "--help" || args[2] == "-h" || args[2] == "-help"))
        throw HelpRequested{ std::string{ DaemonHelp } };
    if (args.size() != 3 || (args[2] != "status" && args[2] != "stop"))
        throw CliError("usage: lima daemon status|stop");
    return Daemon::RunDaemonCommand(ctx, args);
}

namespace {

struct CommandHandlerEntry {
    std::string_view name;
    Cli::CommandHandler handler;
};

constexpr std::array CommandHandlers{
    CommandHandlerEntry{ "mdrun", Cli::RunMdrun },
    CommandHandlerEntry{ "buildmembrane", Cli::RunBuildMembrane },
    CommandHandlerEntry{ "makesimparams", Cli::RunMakeSimParams },
    CommandHandlerEntry{ "selftest", Cli::RunSelfTest },
    CommandHandlerEntry{ "render", Cli::RunRender },
    CommandHandlerEntry{ "makebox", Cli::RunMakeBox },
    CommandHandlerEntry{ "solvate", Cli::RunSolvate },
    CommandHandlerEntry{ "insertmolecule", Cli::RunInsertMolecule },
    CommandHandlerEntry{ "insertmolecules", Cli::RunInsertMolecules },
    CommandHandlerEntry{ "editconf", Cli::RunEditConf },
    CommandHandlerEntry{ "em", Cli::RunEnergyMinimization },
    CommandHandlerEntry{ "togmx", Cli::RunToGmx },
    CommandHandlerEntry{ "daemon", Cli::RunDaemon },
    CommandHandlerEntry{ "wait", Cli::RunWait },
};
static_assert(CommandHandlers.size() == Cli::Commands.size() && std::ranges::all_of(Cli::Commands, [](const auto& command) {
        return std::ranges::find(CommandHandlers, command.name, &CommandHandlerEntry::name) != CommandHandlers.end();
    }), "Every command in Cli::Commands needs exactly one handler");

const Cli::CommandDefinition* FindCommand(const std::string_view name) {
    const auto command = std::ranges::find(Cli::Commands, name, &Cli::CommandDefinition::name);
    return command == Cli::Commands.end() ? nullptr : &*command;
}

Cli::CommandHandler HandlerOf(const Cli::CommandDefinition& command) {
    return std::ranges::find(CommandHandlers, command.name, &CommandHandlerEntry::name)->handler;
}

void PrintGeneralHelp(std::ostream& out) {
    out << "Usage: lima [--no-daemon] COMMAND [OPTION]...\n"
           "       lima [--help] [--version]\n\n"
           "LIMA is a suite of molecular-dynamics and membrane-simulation tools.\n\n"
           "Commands:\n";
    for (const auto& command : Cli::Commands)
        out << Lima::Format("  {:<20} {}\n", command.name, command.summary);
    out << "\nRun 'lima help COMMAND' or 'lima COMMAND --help' for command help.\n";
}

int Dispatch(const Cli::Context& ctx, int argc, char** argv) {
    if (argc == 1) {
        PrintGeneralHelp(ctx.out);
        return 0;
    }

    const std::string_view argument{ argv[1] };
    if (argument == "--help" || argument == "-h" || argument == "-help") {
        PrintGeneralHelp(ctx.out);
        return 0;
    }
    if (argument == "--version") {
        ctx.out << "lima " LIMA_VERSION "\n";
        return 0;
    }
    if (argument == "help") {
        if (argc == 2) {
            PrintGeneralHelp(ctx.out);
            return 0;
        }
        if (argc != 3) throw CliError("usage: lima help COMMAND");
        const auto* command = FindCommand(argv[2]);
        if (!command) throw CliError(Lima::Format("unrecognized command '{}'", argv[2]));
        ctx.out << command->helpText;
        return 0;
    }
    if (argument == "setregistry") {
        if (argc != 2) throw CliError("'setregistry' does not accept options");
        RegisterGrofileAssociation();
        return 0;
    }

    const auto* command = FindCommand(argument);
    if (!command) throw CliError(Lima::Format("unrecognized command '{}'", argument));
    return HandlerOf(*command)(ctx, argc, argv);
}

// Runs a command line in this thread, either in the user's lima process or on a daemon thread
int RunCommand(const Cli::Context& ctx, const std::vector<std::string>& args) {
    std::vector<std::string> storage = args;
    std::vector<char*> argv;
    for (auto& arg : storage) argv.push_back(arg.data());
    try {
        return Dispatch(ctx, static_cast<int>(argv.size()), argv.data());
    }
    catch (const HelpRequested& help) {
        ctx.out << help.helpText;
        return 0;
    }
    catch (const CliError& error) {
        ctx.err << "lima: " << error.what() << "\nTry 'lima --help' for more information.\n";
        return 2;
    }
    catch (const std::exception& error) {
        ctx.err << "lima: " << error.what() << '\n';
        return 1;
    }
    catch (...) {
        ctx.err << "lima: unknown error\n";
        return 1;
    }
}

} // namespace

int main(int argc, char** argv) {
    std::vector<std::string> args(argv, argv + argc);
    if (args.size() >= 2 && args[1] == Cli::Daemon::internalServeCommand) {
        try {
            return Cli::Daemon::Serve(args, RunCommand);
        }
        catch (const std::exception& error) {
            std::cerr << "lima daemon: " << error.what() << '\n';
            return 1;
        }
    }

    const bool noDaemon = args.size() >= 2 && args[1] == "--no-daemon";
    if (noDaemon) args.erase(args.begin() + 1);

    // --async is not an option of mdrun, em or buildmembrane; their handlers never see it. It tells main how to run
    // the command: submit it to the daemon and return at once, instead of waiting for it to finish
    if (const auto async = std::ranges::find(args.begin() + std::min<size_t>(args.size(), 2), args.end(), "--async"); async != args.end()) {
        args.erase(async);
        if (noDaemon || !Cli::Daemon::ShouldUse(args)) {
            std::cerr << "lima: --async runs the command on the lima daemon: only mdrun, em and buildmembrane, without -d,\n"
                         "and not with --no-daemon or LIMA_NO_DAEMON\nTry 'lima --help' for more information.\n";
            return 2;
        }
        try {
            if (const auto exitCode = Cli::Daemon::RunRemote(args, true)) return *exitCode;
        }
        catch (const std::exception& error) {
            std::cerr << "lima: the lima daemon is unavailable (" << error.what() << ")\n";
        }
        std::cerr << "lima: --async needs the lima daemon\n";
        return 1;
    }

    if (!noDaemon && Cli::Daemon::ShouldUse(args)) {
        try {
            if (const auto exitCode = Cli::Daemon::RunRemote(args)) return *exitCode;
        }
        catch (const std::exception& error) {
            std::cerr << "lima: the lima daemon is unavailable (" << error.what() << "), running in this process instead\n";
        }
    }
    return RunCommand(Cli::Context{ fs::current_path(), std::cout, std::cerr }, args);
}
