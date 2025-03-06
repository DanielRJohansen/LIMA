#include "SimParams.h"

#include "Filehandling.h"

#include <unordered_map>
#include <type_traits> // For std::is_integral, std::is_floating_point, and static_assert
#include <functional>

using std::string;
namespace fs = std::filesystem;

using Dictionary = std::unordered_map<string, string>;

void ParseMdp(const Dictionary &mdp_dict, SimParams &params);

template <typename T>
constexpr T convertStringvalueToValue(const std::vector<std::pair<string, T>> pairs, const string& key_str, const string& val_str) {
	for (auto& pair : pairs) {
		if (pair.first == val_str) {
			return pair.second;
		}
	}

	throw std::runtime_error("Illegal key-value pair in sim_params.txt: " + key_str + " " + val_str);
}
// Helper function for overwriting templates
template <typename T>
constexpr void overwriteParamNonNumbers(Dictionary& dict, const std::string& key, T& val, std::function<T(const string&)> transform) {
	if (dict.count(key)) {
		val = transform(dict[key]);
	}
}

void Readb(const Dictionary& dict, bool& value, const std::string& key_name) {
	if (!dict.contains(key_name))
		return;

	std::string value_str = dict.at(key_name);
	if (value_str != "true" && value_str != "false") {
		throw std::runtime_error("Illegal key-value pair in sim_params.txt: " + key_name + " " + value_str);
	}
	
	value = (value_str == "true");
}
template <std::integral T>
constexpr void Readi(const Dictionary& dict, T& param,
	const std::string& key, std::function<T(const T&)> transform = [](const T& v) { return v; })
{
	if (dict.count(key)) {
		param = static_cast<T>(transform(std::stoll(dict.at(key))));
	}
}
template <std::floating_point T, typename Transform = std::function<T(const T&)>>
constexpr void Readf(const Dictionary& dict, T& param,
	const std::string& key, Transform transform = [](const T& v) { return v; })
{
	if (dict.count(key)) {
		param = static_cast<T>(transform(std::stod(dict.at(key))));
	}
}

SimParams::SimParams(const fs::path& path) {
    const bool forceKeysAndValuesLowercase = true;
    auto dict = FileUtils::parseINIFile(path.string(), forceKeysAndValuesLowercase);

    if (path.extension() == ".mdp") {
        ParseMdp(dict, *this);
        return;
    }


    // Main parameters
    Readi(dict, n_steps, "n_steps");
    Readf(dict, dt, "dt", [](auto val) { return val * FEMTO_TO_NANO; }); // [ns]
    Readb(dict, em_variant, "em");
    Readf(dict, em_force_tolerance, "em_force_tolerance");              // [kJ/mol/nm]
    Readi(dict, stepsPerNlistupdate, "stepsPerNlistupdate");

    // Physics parameters
    overwriteParamNonNumbers<BoundaryConditionSelect>(dict, "boundarycondition", bc_select,
        [](const string& value) {
            return convertStringvalueToValue<BoundaryConditionSelect>(
                { {"pbc", PBC}, {"nobc", NoBC} }, "boundarycondition", value);
        }
    );
    Readb(dict, enable_electrostatics, "enable_electrostatics");
    Readf(dict, cutoff_nm, "cutoff_nm"); // [nm]
    // TODO: Parse SupernaturalForcesSelect (SNF) when needed

    // Thermostat parameters
    Readi(dict, steps_per_temperature_measurement, "steps_per_temperature_measurement");
    Readb(dict, apply_thermostat, "apply_thermostat");
    // New thermostat parameters (commented out)
    // Readf(dict, ref_t, "ref_t");         // Reference temperature [K] (critical)
    // Readf(dict, tau_t, "tau_t");         // Temperature coupling constant [ps] (critical)
    // if (dict.count("tcoupl")) {          // Temperature coupling algorithm (important)
    //     tcoupl = dict["tcoupl"];
    // }

    // Barostat (Pressure Coupling) parameters (commented out)
    // if (dict.count("pcoupl")) {          // Pressure coupling algorithm (critical)
    //     pcoupl = dict["pcoupl"];
    // }
    // Readf(dict, ref_p, "ref_p");         // Reference pressure [bar] (critical)
    // Readf(dict, tau_p, "tau_p");         // Pressure coupling constant [ps] (critical)
    // Readf(dict, compressibility, "compressibility"); // Isothermal compressibility [bar^-1] (important)

    // Integration parameters (commented out)
    // if (dict.count("integrator")) {      // Integration algorithm (e.g., md, sd) (critical)
    //     integrator = dict["integrator"];
    // }

    // Constraints (commented out)
    // Readi(dict, constraints, "constraints"); // Constraint type (0: none, 1: bonds, etc.) (important)

    // Electrostatics / PME parameters (commented out)
    // if (dict.count("coulombtype")) {     // Electrostatics method (PME, Reaction Field, etc.) (critical)
    //     coulombtype = dict["coulombtype"];
    // }
    // Readi(dict, pme_order, "pme_order"); // PME interpolation order (important)
    // Readf(dict, fourierspacing, "fourierspacing"); // Fourier grid spacing for PME [nm] (important)

    // Output parameters
    Readi(dict, data_logging_interval, "data_logging_interval");
    Readb(dict, save_energy, "save_energy");
    // TODO: Parse ColoringMethod when needed (e.g., via a string-to-enum conversion)

    // Additional output parameters (commented out)
    // Readi(dict, nstxout, "nstxout");   // Frequency for writing coordinates [steps] (important)
    // Readi(dict, nstvout, "nstvout");   // Frequency for writing velocities [steps] (unimportant)
    // Readi(dict, nstenergy, "nstenergy"); // Frequency for writing energies [steps] (critical)
    // Readi(dict, nstlog, "nstlog");    // Frequency for writing log info [steps] (important)
    // Readi(dict, gen_seed, "gen_seed"); // Random seed for stochastic thermostats (unimportant)
    // Readi(dict, nstcomm, "nstcomm");   // Frequency for center-of-mass removal [steps] (important)
    // if (dict.count("comm_mode")) {     // Center-of-mass removal mode (unimportant)
    //     comm_mode = dict["comm_mode"];
    // }

    // Debug parameters
    Readb(dict, stepwise, "stepwise");
}

void SimParams::dumpToFile(const fs::path& filename) {
    std::ostringstream buffer;

    // Main parameters
    buffer << "\n// Main params\n";
    buffer << "n_steps=" << n_steps << "\n";
    buffer << "dt=" << static_cast<int>(std::round(dt * NANO_TO_FEMTO)) << " # [fs]\n";
    buffer << "em=" << (em_variant ? "true" : "false") << " # Is an energy-minimization sim\n";
    buffer << "em_force_tolerance=" << em_force_tolerance << " # [kJ/mol/nm] - only relevant if em=true\n";
    buffer << "stepsPerNlistupdate=" << stepsPerNlistupdate << " # [steps]\n";

    // Physics parameters
    buffer << "\n// Physics params\n";
    buffer << "boundarycondition=" 
           << (bc_select == PBC ? "PBC" : "NoBC") << " # (PBC, NoBC)\n";
    buffer << "enable_electrostatics=" << (enable_electrostatics ? "true" : "false") << "\n";
    buffer << "cutoff_nm=" << cutoff_nm << " # [nm]\n";
    // buffer << "// SNF parameter not implemented yet\n";

    // Thermostat parameters
    buffer << "\n// Thermostat params\n";
    buffer << "steps_per_temperature_measurement=" << steps_per_temperature_measurement << " # [steps]\n";
    buffer << "apply_thermostat=" << (apply_thermostat ? "true" : "false") 
           << " # Thermostat on/off\n";
    // Uncomment the following lines when needed:
    // buffer << "ref_t=" << ref_t << " # Reference temperature [K] (critical)\n";
    // buffer << "tau_t=" << tau_t << " # Temperature coupling constant [ps] (critical)\n";
    // buffer << "tcoupl=" << tcoupl << " # Temperature coupling algorithm (important)\n";

    // Integration parameters
    buffer << "\n// Integration params\n";
    // Uncomment the following line when needed:
    // buffer << "integrator=" << integrator << " # Integration algorithm (e.g., md, sd) (critical)\n";

    // Constraints
    buffer << "\n// Constraints\n";
    // Uncomment the following line when needed:
    // buffer << "constraints=" << constraints << " # Constraint type (0: none, 1: bonds, etc.) (important)\n";

    // Electrostatics / PME parameters
    buffer << "\n// Electrostatics / PME params\n";
    // Uncomment the following lines when needed:
    // buffer << "coulombtype=" << coulombtype << " # Electrostatics method (PME, Reaction Field, etc.) (critical)\n";
    // buffer << "pme_order=" << pme_order << " # PME interpolation order (important)\n";
    // buffer << "fourierspacing=" << fourierspacing << " # Fourier grid spacing for PME [nm] (important)\n";

    // Output parameters
    buffer << "\n// Output params\n";
    buffer << "data_logging_interval=" << data_logging_interval << " # [steps]\n";
    buffer << "save_energy=" << (save_energy ? "true" : "false") 
           << " # Save kinetic and potential energy to file\n";
    // Uncomment the following lines when needed:
    // buffer << "colormethod=" << (coloring_method == ColoringMethod::Atomname ? "Atomname" : "other")
    //        << " # Coloring method for atoms (important)\n";
    // buffer << "nstxout=" << nstxout << " # Frequency for writing coordinates [steps] (important)\n";
    // buffer << "nstvout=" << nstvout << " # Frequency for writing velocities [steps] (unimportant)\n";
    // buffer << "nstenergy=" << nstenergy << " # Frequency for writing energies [steps] (critical)\n";
    // buffer << "nstlog=" << nstlog << " # Frequency for writing log info [steps] (important)\n";
    // buffer << "gen_seed=" << gen_seed << " # Random seed for stochastic thermostats (unimportant)\n";
    // buffer << "nstcomm=" << nstcomm << " # Frequency for center-of-mass removal [steps] (important)\n";
    // buffer << "comm_mode=" << comm_mode << " # Center-of-mass removal mode (unimportant)\n";

    // Debug parameters
    buffer << "\n// Debug params (for developers)\n";
    buffer << "stepwise=" << (stepwise ? "true" : "false") 
           << " # Wait for user to input key 'N' before each step\n";

    // Write the buffer to file
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Unable to open file: " + filename.string());
    }

    file << buffer.str();
    file.flush();
    file.close();
}



void ParseMdp(const Dictionary &mdp_dict, SimParams &params) {
    // --- Integrator ---
    if (mdp_dict.find("integrator") != mdp_dict.end()) {
        std::string integrator = mdp_dict.at("integrator");
        if (integrator == "md-vv") {
            params.em_variant = false;
        } else if (integrator == "steep") {
            params.em_variant = true;
        } else {
            throw std::runtime_error("Unsupported integrator: " + integrator);
        }
    }
    
    // --- Time & Steps ---
    if (mdp_dict.find("nsteps") != mdp_dict.end()) {
        try {
            params.n_steps = std::stoull(mdp_dict.at("nsteps"));
        } catch (...) {
            throw std::runtime_error("Invalid value for nsteps: " + mdp_dict.at("nsteps"));
        }
    }
    if (mdp_dict.find("dt") != mdp_dict.end()) {
        try {
            double dt_ps = std::stod(mdp_dict.at("dt"));
            params.dt = dt_ps * PICO_TO_NANO; // dt stored in ns
        } catch (...) {
            throw std::runtime_error("Invalid value for dt: " + mdp_dict.at("dt"));
        }
    }
    
    // --- Energy Minimization ---
    if (mdp_dict.find("emtol") != mdp_dict.end()) {
        try {
            params.em_force_tolerance = std::stod(mdp_dict.at("emtol"));
        } catch (...) {
            throw std::runtime_error("Invalid value for emtol: " + mdp_dict.at("emtol"));
        }
    }
    if (mdp_dict.find("nstlist") != mdp_dict.end()) {
        try {
            params.stepsPerNlistupdate = std::stoi(mdp_dict.at("nstlist"));
        } catch (...) {
            throw std::runtime_error("Invalid value for nstlist: " + mdp_dict.at("nstlist"));
        }
    }
    
    // --- Physics ---
    // Boundary condition: using mdp key "pbc"
    if (mdp_dict.find("pbc") != mdp_dict.end()) {
        std::string pbc = mdp_dict.at("pbc");
        if (pbc == "xyz") {
            params.bc_select = PBC;
        } else if (pbc == "none" || pbc == "no") {
            params.bc_select = NoBC;
        } else {
            throw std::runtime_error("Unsupported pbc value: " + pbc);
        }
    }
    
    // Electrostatics: only support PME
    if (mdp_dict.find("coulombtype") != mdp_dict.end()) {
        std::string coulombtype = mdp_dict.at("coulombtype");
        if (coulombtype != "PME") {
            throw std::runtime_error("Unsupported coulombtype: " + coulombtype);
        }
        // Set enable_electrostatics true (default is true).
        params.enable_electrostatics = true;
    }
    
    // Cutoff: using mdp key "rcoulomb"
    if (mdp_dict.find("rcoulomb") != mdp_dict.end()) {
        try {
            params.cutoff_nm = std::stof(mdp_dict.at("rcoulomb"));
        } catch (...) {
            throw std::runtime_error("Invalid value for rcoulomb: " + mdp_dict.at("rcoulomb"));
        }
    }
    
    // --- Output Intervals ---
    std::vector<std::string> outKeys = { "nstxout", "nstvout", "nstenergy", "nstlog", "nstcomm" };
    int commonInterval = -1;
    bool found = false;
    for (const auto &key : outKeys) {
        if (mdp_dict.find(key) != mdp_dict.end()) {
            int interval;
            try {
                interval = std::stoi(mdp_dict.at(key));
            } catch (...) {
                throw std::runtime_error("Invalid integer for " + key + ": " + mdp_dict.at(key));
            }
            if (!found) {
                commonInterval = interval;
                found = true;
            } else if (interval != commonInterval) {
                throw std::runtime_error("Output interval mismatch for key '" + key +
                                         "': value " + std::to_string(interval) +
                                         " does not match expected " + std::to_string(commonInterval));
            }
        }
    }
    if (found) {
        params.data_logging_interval = commonInterval;
    }
    
    // --- Thermostat ---
    // Determine whether thermostat is applied from mdp key "tcoupl".
    if (mdp_dict.find("tcoupl") != mdp_dict.end()) {
        std::string tcoupl = mdp_dict.at("tcoupl");
        if (tcoupl == "no" || tcoupl == "none") {
            params.apply_thermostat = false;
        } else {
            params.apply_thermostat = true;
        }
    }
    // Frequency of temperature coupling: "nsttcouple"
    if (mdp_dict.find("nsttcouple") != mdp_dict.end()) {
        try {
            params.steps_per_temperature_measurement = std::stoll(mdp_dict.at("nsttcouple"));
        } catch (...) {
            throw std::runtime_error("Invalid value for nsttcouple: " + mdp_dict.at("nsttcouple"));
        }
    }
}
