#include <iostream>
#include <format>
#include <filesystem>

#include <CommandlineUtils.h>

#include "Environment.h"
namespace fs = std::filesystem;

struct MdrunSetup {
    MdrunSetup(int argc, char** argv) : work_dir(fs::current_path()) {
        for (int i = 2; i < argc; ++i) {
            std::string arg = CmdLineUtils::ToLowercase(argv[i]);

            if (arg == "-conf") {
                if (i + 1 < argc) {
                    conf = work_dir / argv[++i];
                }
                else {
                    std::cerr << "-conf expects a path argument." << std::endl;
                    exit(1);
                }
            }
            else if (arg == "-topology") {
                if (i + 1 < argc) {
                    topol = work_dir / argv[++i];
                }
                else {
                    std::cerr << "-topology expects a path argument." << std::endl;
                    exit(1);
                }
            }
            else if (arg == "-simparams") {
                if (i + 1 < argc) {
                    simpar = work_dir / argv[++i];
                }
                else {
                    std::cerr << "-simparams expects a path argument." << std::endl;
                    exit(1);
                }
            }
            else if (arg == "-structure_out") {
                if (i + 1 < argc) {
                    conf_out = work_dir / argv[++i];
                }
                else {
                    std::cerr << "-structure_out expects a path argument." << std::endl;
                    exit(1);
                }
            }
            else if (arg == "-help" || arg == "-h") {
                std::cout << helpText;
                exit(0);
            }
            else if (arg == "-display" || arg == "-d")
                envmode = Full;
            else {
                std::cerr << "Unknown argument: " << arg << std::endl;
                exit(1);
            }
        }
    }

    EnvMode envmode = ConsoleOnly;
    fs::path work_dir;
    fs::path conf = work_dir / "molecule/conf.gro";
    fs::path topol = work_dir / "molecule/topol.top";
    fs::path simpar = work_dir / "sim_params.txt";
    fs::path conf_out = work_dir / "out.gro";

private:
    const std::string helpText = R"(
Usage: mdrun [OPTIONS]

Description:
    This program runs a molecular dynamics simulation based on provided configuration files and parameters.

Options:
    -conf [path]
        Path to the configuration (.gro) file. Defaults to molecule/conf.gro.
    
    -topology [path]
        Path to the topology (.top) file. Defaults to molecule/topol.top.
    
    -simparams [path]
        Path to the simulation parameters file. Defaults to sim_params.txt.

    -display, -d
		Flag to enable the display, rendering the simulation and displaying information such as temperature, step and more.

    -structure_out [path]
        Output path for the resulting structure file. Defaults to out.gro.

    -help, -h
        Display this help text and exit.

Example:
    mdrun -conf myconf.gro -topology mytopol.top -simparams params.txt -display
    )";
};

int mdrun(int argc, char** argv) {
    auto t0 = std::chrono::steady_clock::now();
    MdrunSetup setup(argc, argv);
    auto env = std::make_unique<Environment>(setup.work_dir, setup.envmode);

    const SimParams ip(setup.simpar);
    GroFile grofile{ setup.conf };
    TopologyFile topfile{ setup.topol };

    env->CreateSimulation(grofile, topfile, ip);
    env->run();

    env->WriteBoxCoordinatesToFile(grofile);
    grofile.printToFile(setup.conf_out);


    // TODO: Move this to a helper namespace somewhere

    // Calculate total time simulated (in nanoseconds)
    const double total_ns = static_cast<double>(ip.n_steps) * ip.dt * LIMA_TO_NANO;

    // Measure elapsed time in seconds
    auto t1 = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = t1 - t0;
    double elapsed_time_sec = elapsed_seconds.count();

    // Calculate performance metrics
    double ns_per_s = total_ns / elapsed_time_sec;
    double ns_per_hr = ns_per_s * 3600.0;

    // Print time and performance info
    printf("Total simulation time: %.2f seconds\n", elapsed_time_sec);
    printf("Performance: %.3f ns/s, %.3f ns/hr\n", ns_per_s, ns_per_hr);


    return 0;
}
