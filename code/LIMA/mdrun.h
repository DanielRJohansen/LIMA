#include <iostream>
#include <format>
#include <filesystem>
#include <CommandlineUtils.h>

#include "Environment.h"
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

int mdrun(int argc, char** argv) {

    const std::string helpText = R"(
Usage: mdrun [OPTIONS]

Description:
    This program runs a molecular dynamics simulation based on provided configuration files and parameters.

Options:
    -conf [path]
        Path to the configuration (.gro) file. Defaults to ./conf.gro
    
    -topology [path]
        Path to the topology (.top) file. Defaults to ./topol.top
    
    -simparams [path]
        Path to the simulation parameters file. Defaults to ./sim_params.txt

    -display, -d
        Flag to enable the display, rendering the simulation and displaying information such as temperature, step and more

    -conf_out [path]
        Output path for the resulting structure file. Defaults to ./out.gro

    -trajectory, -trr [path]
		Output path for the trajectory file. Defaults to not outputting a trajectory file

    -help, -h
        Display this help text and exit.

Example:
    mdrun -conf myconf.gro -topology mytopol.top -simparams params.txt -display
    )";

    ArgParser parser(helpText);

    fs::path work_dir;
    fs::path conf = work_dir / "conf.gro";
    fs::path topol = work_dir / "topol.top";
    fs::path simpar = work_dir / "sim_params.txt";
    fs::path conf_out = work_dir / "out.gro";
    fs::path trajOut{};

    bool render = false;

    parser.AddOption({ "-conf", "-c" }, false, conf);
	parser.AddOption({ "-topology", "-top", "-t"}, false, topol);
	parser.AddOption({ "-simparams", "-s" }, false, simpar);
	parser.AddOption({ "-conf_out", "-co" }, false, conf_out);
	parser.AddOption({ "-trajectory", "-trr", "-traj"}, false, trajOut);
	parser.AddFlag({ "-display", "-d" }, [&render]() { render = true; });

    EnvMode envmode = render ? Full : ConsoleOnly;

    auto env = std::make_unique<Environment>(work_dir, envmode);

    const SimParams ip(simpar);
    GroFile grofile{ conf };
    TopologyFile topfile{ topol };

    if ((NodeIndex(grofile.box_size.ToInt3()).toFloat3() - grofile.box_size).len() > 0.0001 || grofile.box_size.x != grofile.box_size.y || grofile.box_size.y != grofile.box_size.z) {
        const int newSize = std::ceil(std::max(std::max(grofile.box_size.x, grofile.box_size.y), grofile.box_size.z));
        printf("Boxsize was not an integer, or was not cubic. Setting new boxsize to %d", newSize);
        grofile.box_size = Float3(newSize, newSize, newSize);
    }

    env->CreateSimulation(grofile, topfile, ip);

    const std::chrono::duration<double> enginetime = env->run(); // [s]
    printf("Engine time %f\n", enginetime.count());
    env->WriteBoxCoordinatesToFile(grofile);
    grofile.printToFile(conf_out);

    if (!trajOut.empty()) {
        env->getSimPtr()->ToTracjectoryFile()->Dump(trajOut);
    }

    // Calculate total time simulated (in nanoseconds)
    const double total_ns_simulated = static_cast<double>(ip.n_steps) * ip.dt;
	PrintTiming(enginetime, total_ns_simulated);

    return 0;
}
