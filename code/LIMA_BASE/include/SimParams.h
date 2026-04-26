#pragma once

#include "LimaTypes.cuh"
#include <set>
#include <filesystem>

enum class ColoringMethod { Atomname, Charge, GradientFromAtomid, PersistentClusterId, ForceMagnitude };

enum BoundaryConditionSelect{NoBC, PBC};

enum SupernaturalForcesSelect{None, HorizontalSqueeze, HorizontalChargeField, BoxEdgePotential, ElasticPosition};

struct SimParams {
    SimParams() {}
    SimParams(const std::filesystem::path& path);
    SimParams(std::initializer_list<int>) = delete;

    void dumpToFile(const std::filesystem::path& filename = "sim_params.txt");

    // Main parameters
    uint64_t n_steps = 1000;
    float dt = 2.f * FEMTO_TO_NANO;           // Time step [ns]
    bool em_variant = false;
    float em_force_tolerance = 1000;           // [kJ/mol/nm]
    int stepsPerNlistupdate = 5;

    // Physics parameters
    BoundaryConditionSelect bc_select{ PBC };
    bool enable_electrostatics = true;
    float cutoff_nm = 1.2f;                    // Cutoff distance [nm]
    std::set<SupernaturalForcesSelect> snf_select;

    // Thermostat
    int64_t steps_per_temperature_measurement = 200;
    bool apply_thermostat = false;
    float ref_t = 300.0f;                      // Reference temperature [K] (critical)
    // float tau_t = 0.1f;                     // Temperature coupling constant [ps] (critical)
    // std::string tcoupl = "V-rescale";         // Temperature coupling algorithm (important)

    // Barostat (Pressure Coupling)
    // std::string pcoupl = "Parrinello-Rahman"; // Pressure coupling algorithm (critical)
    // float ref_p = 1.0f;                     // Reference pressure [bar] (critical)
    // float tau_p = 2.0f;                     // Pressure coupling constant [ps] (critical)
    // float compressibility = 4.5e-5f;        // Isothermal compressibility [bar^-1] (important)

    // Integration parameters
    // std::string integrator = "md";          // Integration algorithm (e.g., md, sd) (critical)

    // Constraints
    // int constraints = 0;                    // Constraint type (0: none, 1: bonds, etc.) (important)

    // Electrostatics / PME parameters
    // std::string coulombtype = "PME";        // Electrostatics method (PME, Reaction Field, etc.) (critical)
    // int pme_order = 4;                      // PME interpolation order (important)
    // float fourierspacing = 0.12f;           // Fourier grid spacing for PME [nm] (important)

    // Output parameters
    int data_logging_interval = 5;
    bool save_energy = false;
    ColoringMethod coloring_method = ColoringMethod::Atomname;  // TODO: THis is actually being ignored now...
    // int nstxout = 500;                    // Frequency for writing coordinates [steps] (important)
    // int nstvout = 500;                    // Frequency for writing velocities [steps] (unimportant)
    // int nstenergy = 100;                  // Frequency for writing energies [steps] (critical)
    // int nstlog = 100;                     // Frequency for writing log info [steps] (important)
    // unsigned int gen_seed = 42;           // Random seed for stochastic thermostats (unimportant)
    // int nstcomm = 10;                     // Frequency for center-of-mass removal [steps] (important)
    // std::string comm_mode = "linear";     // Center-of-mass removal mode (unimportant)

    // Debug parameters
    bool stepwise = false;                   // Wait for user input key "N" before each step
};
