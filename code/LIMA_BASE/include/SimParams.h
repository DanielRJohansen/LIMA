#pragma once

#include "LimaTypes.cuh"

#include <filesystem>
#include <format>
#include <optional>
#include <set>
#include <string_view>
#include <tuple>
#include <utility>

enum class ColoringMethod { Atomname, Charge, GradientFromAtomid, PersistentClusterId, ForceMagnitude, NewCartoon };
enum BoundaryConditionSelect { NoBC, PBC };
enum SupernaturalForcesSelect { None, HorizontalSqueeze, HorizontalChargeField, BoxEdgePotential, ElasticPosition };
enum class SimParamSection { Main, Physics, Thermostat, Output, Debug };

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

template<typename T, typename CharT>
struct std::formatter<SimParam<T>, CharT> : std::formatter<T, CharT> {
    template<typename FormatContext>
    auto format(const SimParam<T>& param, FormatContext& context) const {
        return std::formatter<T, CharT>::format(param.value, context);
    }
};

struct SimParams {
    SimParams() = default;
    explicit SimParams(const std::filesystem::path& path);
    SimParams(std::initializer_list<int>) = delete;

    void DumpToFile(const std::filesystem::path& filename = "sim_params.txt") const;

    SimParam<uint64_t> n_steps{ "n_steps", 1000, SimParamSection::Main };
    SimParam<float> dt{ "dt", 2.f * FEMTO_TO_NANO, SimParamSection::Main, "[fs]" };
    SimParam<bool> em_variant{ "em", false, SimParamSection::Main, "Is an energy-minimization simulation" };
    SimParam<float> em_force_tolerance{ "em_force_tolerance", 1000.f, SimParamSection::Main,
        "[kJ/mol/nm], only relevant if em=true" };
    SimParam<int> stepsPerNlistupdate{ "stepsPerNlistupdate", 20, SimParamSection::Main, "[steps]" };

    SimParam<BoundaryConditionSelect> bc_select{ "boundarycondition", PBC, SimParamSection::Physics, "PBC or NoBC" };
    SimParam<bool> enable_electrostatics{ "enable_electrostatics", true, SimParamSection::Physics };
    SimParam<float> cutoff_nm{ "cutoff_nm", 1.2f, SimParamSection::Physics, "[nm]" };
    SimParam<std::set<SupernaturalForcesSelect>> snf_select{ "supernatural_forces", {}, SimParamSection::Physics };

    SimParam<int64_t> steps_per_temperature_measurement{ "steps_per_temperature_measurement", 200,
        SimParamSection::Thermostat, "[steps]" };
    SimParam<bool> apply_thermostat{ "apply_thermostat", false, SimParamSection::Thermostat, "Thermostat on/off" };
    SimParam<float> ref_t{ "ref_t", 300.f, SimParamSection::Thermostat, "Reference temperature [K]" };

    SimParam<int> data_logging_interval{ "data_logging_interval", 5, SimParamSection::Output, "[steps]" };
    SimParam<bool> save_energy{ "save_energy", false, SimParamSection::Output,
        "Save kinetic and potential energy to file" };
    SimParam<ColoringMethod> coloring_method{ "coloring_method", ColoringMethod::Atomname, SimParamSection::Output };

    SimParam<bool> stepwise{ "stepwise", false, SimParamSection::Debug,
        "Wait for user to input key 'N' before each step" };

    auto Params() {
        return std::tie(n_steps, dt, em_variant, em_force_tolerance, stepsPerNlistupdate,
            bc_select, enable_electrostatics, cutoff_nm, snf_select,
            steps_per_temperature_measurement, apply_thermostat, ref_t,
            data_logging_interval, save_energy, coloring_method, stepwise);
    }

    auto Params() const {
        return std::tie(n_steps, dt, em_variant, em_force_tolerance, stepsPerNlistupdate,
            bc_select, enable_electrostatics, cutoff_nm, snf_select,
            steps_per_temperature_measurement, apply_thermostat, ref_t,
            data_logging_interval, save_energy, coloring_method, stepwise);
    }
};
