#pragma once

#include "LimaTypes.cuh"

#include <filesystem>
#ifndef __CUDACC__
#include <format>
#endif
#include <optional>
#include <set>
#include <string_view>
#include <tuple>
#include <utility>

enum class ColoringMethod { Atomname, Charge, GradientFromAtomid, PersistentClusterId, ForceMagnitude, NewCartoon, Lipids, Licorice, Surface };
enum BoundaryConditionSelect { NoBC, PBC };
enum class ConstraintsSelect { None, HBonds };
enum SupernaturalForcesSelect { None, HorizontalSqueeze, HorizontalChargeField, BoxEdgePotential, ElasticPosition };
enum class SimParamSection { Main, Physics, Thermostat, Barostat, Output, Analysis, Debug };

template<typename T>
struct SimParam {
    using Type = T;

    constexpr SimParam(std::string_view name, T defaultValue, SimParamSection section,
        std::optional<std::string_view> comment = std::nullopt)
        : name(name), defaultValue(defaultValue), value(std::move(defaultValue)), section(section), comment(comment) {}

    constexpr operator const T&() const { return value; }

    constexpr SimParam& operator=(const T& newValue) {
        value = newValue;
        return *this;
    }

    constexpr auto empty() const requires requires(const T& candidate) { candidate.empty(); } {
        return value.empty();
    }

    template<typename... Args>
    constexpr decltype(auto) insert(Args&&... args)
        requires requires(T& candidate) { candidate.insert(std::forward<Args>(args)...); } {
        return value.insert(std::forward<Args>(args)...);
    }

    template<typename... Args>
    constexpr decltype(auto) erase(Args&&... args)
        requires requires(T& candidate) { candidate.erase(std::forward<Args>(args)...); } {
        return value.erase(std::forward<Args>(args)...);
    }

    template<typename... Args>
    constexpr decltype(auto) contains(Args&&... args) const
        requires requires(const T& candidate) { candidate.contains(std::forward<Args>(args)...); } {
        return value.contains(std::forward<Args>(args)...);
    }

    std::string_view name;
    T defaultValue;
    T value;
    SimParamSection section;
    std::optional<std::string_view> comment;
};

// Not for CUDA files, which have no use for it: <format> takes seconds to compile, three times over in nvcc
#ifndef __CUDACC__
template<typename T, typename CharT>
struct std::formatter<SimParam<T>, CharT> : std::formatter<T, CharT> {
    template<typename FormatContext>
    auto format(const SimParam<T>& param, FormatContext& context) const {
        return std::formatter<T, CharT>::format(param.value, context);
    }
};
#endif

struct SimParams {
    SimParams() = default;
    explicit SimParams(const std::filesystem::path& path);
    SimParams(std::initializer_list<int>) = delete;

    // Generate the default SimParams for basic EM.
    static SimParams BasicEMSimParams(float forceTolerance = 100.f);

    void DumpToFile(const std::filesystem::path& filename = "sim_params.txt") const;
    // Sets a parameter as a sim_params file line "key=value" would. Throws on an unknown key
    void Set(std::string_view key, std::string_view value);

    SimParam<uint64_t> n_steps{ "n_steps", 1000, SimParamSection::Main };
    SimParam<float> dt{ "dt", 2.f * FEMTO_TO_NANO, SimParamSection::Main, "[fs]" };
    SimParam<bool> em_variant{ "em", false, SimParamSection::Main, "Is an energy-minimization simulation" };
    SimParam<float> em_force_tolerance{ "em_force_tolerance", 1000.f, SimParamSection::Main,
        "[kJ/mol/nm], only relevant if em=true" };
    SimParam<int> stepsPerNlistupdate{ "stepsPerNlistupdate", 20, SimParamSection::Main, "[steps]" };
    SimParam<float> verlet_buffer_tolerance{ "verlet_buffer_tolerance", -1.f, SimParamSection::Main,
        "[kJ/mol/ps] per atom. When positive, the pair list buffer is set so pairs missed between list updates cause at most this "
        "energy drift, as GROMACS' verlet-buffer-tolerance (eg. 0.005). That needs about twice the default buffer, which is ~11% "
        "slower without dynamic pruning. -1 uses rlist" };
    SimParam<float> rlist{ "rlist", 0.f, SimParamSection::Main,
        "[nm] Pair list radius with verlet_buffer_tolerance=-1, at least the cutoff. 0 is the cutoff + 0.05" };

    SimParam<BoundaryConditionSelect> bc_select{ "boundarycondition", PBC, SimParamSection::Physics, "PBC or NoBC" };
    SimParam<bool> enable_electrostatics{ "enable_electrostatics", true, SimParamSection::Physics };
    SimParam<float> cutoff_nm{ "cutoff_nm", 1.2f, SimParamSection::Physics, "[nm]" };
    SimParam<float> fourier_spacing{ "fourierspacing", 0.f, SimParamSection::Physics,
        "[nm] Largest PME grid spacing. 0 is 0.125 * cutoff_nm, the accuracy of GROMACS' defaults (fourierspacing 0.12 "
        "at a 1 nm cutoff), as the Ewald splitting, and with it the grid needed, scales with the cutoff" };
    SimParam<std::set<SupernaturalForcesSelect>> snf_select{ "supernatural_forces", {}, SimParamSection::Physics };
    SimParam<ConstraintsSelect> constraints{ "constraints", ConstraintsSelect::None, SimParamSection::Physics,
        "none or h-bonds: h-bonds holds bonds to hydrogens at their equilibrium length and makes water rigid, as GROMACS' "
        "constraints = h-bonds with settled water. Only 3-atom water is made rigid, eg. TIP4P only gets its O-H bonds constrained. "
        "MD only, energy minimization keeps them flexible" };

    SimParam<int64_t> steps_per_temperature_measurement{ "steps_per_temperature_measurement", 200,
        SimParamSection::Thermostat, "[steps]" };
    SimParam<bool> apply_thermostat{ "apply_thermostat", false, SimParamSection::Thermostat,
        "v-rescale thermostat (Bussi 2007, as GROMACS' tcoupl = v-rescale) on/off. Couples every steps_per_temperature_measurement steps" };
    SimParam<float> ref_t{ "ref_t", 300.f, SimParamSection::Thermostat, "Reference temperature [K]" };
    SimParam<float> tau_t{ "tau_t", 0.1f, SimParamSection::Thermostat, "[ps] Time constant of the thermostat" };
    SimParam<int64_t> thermostat_seed{ "thermostat_seed", 1, SimParamSection::Thermostat,
        "Seed of the thermostat's random numbers. Runs with the same seed and input are identical" };
    SimParam<float> initial_temperature{ "initial_temperature", 0.f, SimParamSection::Thermostat,
        "[K] Give every atom a random velocity at this temperature before the first step, eg. for the first run after "
        "energy minimization. 0 keeps the velocities of the input coordinates, as when continuing a run" };
    SimParam<int64_t> velocity_seed{ "velocity_seed", 1, SimParamSection::Thermostat,
        "Seed of the initial_temperature velocities, as GROMACS' gen-seed. Runs with the same seed and input are identical" };

    SimParam<bool> apply_barostat{ "apply_barostat", false, SimParamSection::Barostat,
        "Isotropic Berendsen pressure coupling on/off" };
    SimParam<int64_t> steps_per_pressure_measurement{ "steps_per_pressure_measurement", 20,
        SimParamSection::Barostat, "[steps]. Align with stepsPerNlistupdate to reuse its neighbor-list rebuild" };
    SimParam<float> ref_p{ "ref_p", 1.f, SimParamSection::Barostat, "Reference pressure [bar]" };
    SimParam<float> tau_p{ "tau_p", 2.f, SimParamSection::Barostat, "[ps] Time constant of the barostat" };
    SimParam<float> compressibility{ "compressibility", 4.5e-5f, SimParamSection::Barostat, "[bar^-1]" };

    SimParam<int> data_logging_interval{ "data_logging_interval", 5, SimParamSection::Output, "[steps]" };
    SimParam<bool> save_energy{ "save_energy", false, SimParamSection::Output,
        "Write energy.xvg: potential, kinetic and total energy [kJ/mol] and temperature [K] every data_logging_interval steps, "
        "in the format of GROMACS' gmx energy" };
    SimParam<ColoringMethod> coloring_method{ "coloring_method", ColoringMethod::Atomname, SimParamSection::Output };

    SimParam<bool> self_diffusion{ "self_diffusion", false, SimParamSection::Analysis,
        "Self-diffusion of the water oxygens, written to self_diffusion.csv. Skips the first 20% of the run" };
    SimParam<bool> density_profile{ "density_profile", false, SimParamSection::Analysis,
        "Water, lipid head and tail density along z, written to density_profile.csv" };

    SimParam<bool> stepwise{ "stepwise", false, SimParamSection::Debug,
        "Wait for user to input key 'N' before each step" };

    auto Params() {
        return std::tie(n_steps, dt, em_variant, em_force_tolerance, stepsPerNlistupdate, verlet_buffer_tolerance, rlist,
            bc_select, enable_electrostatics, cutoff_nm, fourier_spacing, snf_select, constraints,
            steps_per_temperature_measurement, apply_thermostat, ref_t, tau_t, thermostat_seed, initial_temperature, velocity_seed,
            apply_barostat, steps_per_pressure_measurement, ref_p, tau_p, compressibility,
            data_logging_interval, save_energy, coloring_method, self_diffusion, density_profile, stepwise);
    }

    auto Params() const {
        return std::tie(n_steps, dt, em_variant, em_force_tolerance, stepsPerNlistupdate, verlet_buffer_tolerance, rlist,
            bc_select, enable_electrostatics, cutoff_nm, fourier_spacing, snf_select, constraints,
            steps_per_temperature_measurement, apply_thermostat, ref_t, tau_t, thermostat_seed, initial_temperature, velocity_seed,
            apply_barostat, steps_per_pressure_measurement, ref_p, tau_p, compressibility,
            data_logging_interval, save_energy, coloring_method, self_diffusion, density_profile, stepwise);
    }
};
