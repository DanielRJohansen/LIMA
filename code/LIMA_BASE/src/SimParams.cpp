#include "SimParams.h"

#include "Filehandling.h"

#include <algorithm>
#include <array>
#include <charconv>
#include <cctype>
#include <concepts>
#include <format>
#include <fstream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>

namespace fs = std::filesystem;
using Dictionary = std::unordered_map<std::string, std::string>;

namespace {

std::string Lowercase(std::string_view value) {
    std::string result{ value };
    std::ranges::transform(result, result.begin(), [](const unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return result;
}

template<typename T>
T ParseNumber(const std::string_view key, const std::string_view input) {
    T value{};
    const char* begin = input.data();
    const char* end = begin + input.size();
    const auto result = std::from_chars(begin, end, value);
    if (result.ec != std::errc{} || result.ptr != end)
        throw std::runtime_error(std::format("Invalid value for '{}': '{}'", key, input));
    return value;
}

template<typename T>
void ParseValue(const std::string_view key, const std::string_view input, T& value) {
    if constexpr (std::same_as<T, bool>) {
        if (input == "true") value = true;
        else if (input == "false") value = false;
        else throw std::runtime_error(std::format("Invalid boolean for '{}': '{}'", key, input));
    }
    else if constexpr (std::integral<T> || std::floating_point<T>) {
        value = ParseNumber<T>(key, input);
    }
    else if constexpr (std::same_as<T, BoundaryConditionSelect>) {
        if (input == "pbc") value = PBC;
        else if (input == "nobc") value = NoBC;
        else throw std::runtime_error(std::format("Invalid boundary condition: '{}'", input));
    }
    else if constexpr (std::same_as<T, ColoringMethod>) {
        static const std::unordered_map<std::string_view, ColoringMethod> values{
            { "atomname", ColoringMethod::Atomname },
            { "charge", ColoringMethod::Charge },
            { "gradientfromatomid", ColoringMethod::GradientFromAtomid },
            { "persistentclusterid", ColoringMethod::PersistentClusterId },
            { "forcemagnitude", ColoringMethod::ForceMagnitude },
            { "newcartoon", ColoringMethod::NewCartoon }
        };
        const auto found = values.find(input);
        if (found == values.end()) throw std::runtime_error(std::format("Invalid coloring method: '{}'", input));
        value = found->second;
    }
    else if constexpr (std::same_as<T, std::set<SupernaturalForcesSelect>>) {
        value.clear();
        std::istringstream stream{ std::string{ input } };
        for (std::string item; std::getline(stream, item, ',');) {
            if (item.empty() || item == "none") continue;
            if (item == "horizontalsqueeze") value.insert(HorizontalSqueeze);
            else if (item == "horizontalchargefield") value.insert(HorizontalChargeField);
            else if (item == "boxedgepotential") value.insert(BoxEdgePotential);
            else if (item == "elasticposition") value.insert(ElasticPosition);
            else throw std::runtime_error(std::format("Invalid supernatural force: '{}'", item));
        }
    }
}

std::string FormatValue(const bool value) { return value ? "true" : "false"; }
std::string FormatValue(const BoundaryConditionSelect value) { return value == PBC ? "PBC" : "NoBC"; }

std::string FormatValue(const ColoringMethod value) {
    switch (value) {
    case ColoringMethod::Atomname: return "Atomname";
    case ColoringMethod::Charge: return "Charge";
    case ColoringMethod::GradientFromAtomid: return "GradientFromAtomid";
    case ColoringMethod::PersistentClusterId: return "PersistentClusterId";
    case ColoringMethod::ForceMagnitude: return "ForceMagnitude";
    case ColoringMethod::NewCartoon: return "NewCartoon";
    }
    throw std::runtime_error("Invalid coloring method");
}

std::string FormatValue(const std::set<SupernaturalForcesSelect>& values) {
    if (values.empty()) return "None";
    std::string result;
    const auto append = [&result](const std::string_view value) {
        if (!result.empty()) result += ',';
        result += value;
    };
    for (const auto value : values) {
        switch (value) {
        case None: break;
        case HorizontalSqueeze: append("HorizontalSqueeze"); break;
        case HorizontalChargeField: append("HorizontalChargeField"); break;
        case BoxEdgePotential: append("BoxEdgePotential"); break;
        case ElasticPosition: append("ElasticPosition"); break;
        }
    }
    return result.empty() ? "None" : result;
}

std::string_view SectionName(const SimParamSection section) {
    switch (section) {
    case SimParamSection::Main: return "Main parameters";
    case SimParamSection::Physics: return "Physics parameters";
    case SimParamSection::Thermostat: return "Thermostat parameters";
    case SimParamSection::Output: return "Output parameters";
    case SimParamSection::Debug: return "Debug parameters";
    }
    throw std::runtime_error("Invalid simulation parameter section");
}

template<typename T>
std::string FormatValue(const T value) {
    return std::format("{}", value);
}

template<typename Function>
void ForEachParam(SimParams& params, Function&& function) {
    std::apply([&](auto&... param) { (function(param), ...); }, params.Params());
}

template<typename Function>
void ForEachParam(const SimParams& params, Function&& function) {
    std::apply([&](const auto&... param) { (function(param), ...); }, params.Params());
}

void ParseMdp(const Dictionary& mdp, SimParams& params) {
    const auto parseIfPresent = [&mdp]<typename T>(const std::string_view key, T& target) {
        if (const auto found = mdp.find(std::string{ key }); found != mdp.end())
            ParseValue(key, found->second, target);
    };

    if (const auto found = mdp.find("integrator"); found != mdp.end()) {
        if (found->second == "md-vv") params.em_variant.value = false;
        else if (found->second == "steep") params.em_variant.value = true;
        else throw std::runtime_error("Unsupported integrator: " + found->second);
    }

    parseIfPresent("nsteps", params.n_steps.value);
    if (const auto found = mdp.find("dt"); found != mdp.end())
        params.dt.value = ParseNumber<float>("dt", found->second) * PICO_TO_NANO;
    parseIfPresent("emtol", params.em_force_tolerance.value);
    parseIfPresent("nstlist", params.stepsPerNlistupdate.value);

    if (const auto found = mdp.find("pbc"); found != mdp.end()) {
        if (found->second == "xyz") params.bc_select.value = PBC;
        else if (found->second == "none" || found->second == "no") params.bc_select.value = NoBC;
        else throw std::runtime_error("Unsupported pbc value: " + found->second);
    }
    if (const auto found = mdp.find("coulombtype"); found != mdp.end()) {
        if (found->second != "pme") throw std::runtime_error("Unsupported coulombtype: " + found->second);
        params.enable_electrostatics.value = true;
    }
    parseIfPresent("rcoulomb", params.cutoff_nm.value);

    constexpr std::array outputKeys{ "nstxout", "nstvout", "nstenergy", "nstlog", "nstcomm" };
    std::optional<int> outputInterval;
    for (const std::string_view key : outputKeys) {
        if (const auto found = mdp.find(std::string{ key }); found != mdp.end()) {
            const int interval = ParseNumber<int>(key, found->second);
            if (outputInterval && *outputInterval != interval)
                throw std::runtime_error(std::format("Output interval '{}' differs from earlier intervals", key));
            outputInterval = interval;
        }
    }
    if (outputInterval) params.data_logging_interval.value = *outputInterval;

    if (const auto found = mdp.find("tcoupl"); found != mdp.end())
        params.apply_thermostat.value = found->second != "no" && found->second != "none";
    parseIfPresent("nsttcouple", params.steps_per_temperature_measurement.value);
}

} // namespace

SimParams::SimParams(const fs::path& path) {
    const Dictionary dictionary = FileUtils::parseINIFile(path.string(), true);
    if (path.extension() == ".mdp") {
        ParseMdp(dictionary, *this);
        return;
    }

    std::unordered_set<std::string> knownKeys;
    ForEachParam(*this, [&](auto& param) {
        const std::string key = Lowercase(param.name);
        knownKeys.insert(key);
        if (const auto found = dictionary.find(key); found != dictionary.end()) {
            ParseValue(param.name, found->second, param.value);
            if constexpr (std::same_as<typename std::remove_cvref_t<decltype(param)>::Type, float>)
                if (param.name == "dt") param.value *= FEMTO_TO_NANO;
        }
    });

    for (const auto& [key, value] : dictionary)
        if (!knownKeys.contains(key))
            throw std::runtime_error(std::format("Unknown simulation parameter '{}={}'", key, value));
}

void SimParams::DumpToFile(const fs::path& filename) const {
    std::ofstream file{ filename };
    if (!file) throw std::runtime_error("Unable to open file: " + filename.string());

    std::optional<SimParamSection> currentSection;
    ForEachParam(*this, [&](const auto& param) {
        if (currentSection != param.section) {
            if (currentSection) file << '\n';
            currentSection = param.section;
            file << "// " << SectionName(param.section) << '\n';
        }

        file << param.name << '=';
        if constexpr (std::same_as<typename std::remove_cvref_t<decltype(param)>::Type, float>) {
            if (param.name == "dt") {
                file << FormatValue(static_cast<float>(param.value * NANO_TO_FEMTO));
            }
            else {
                file << FormatValue(param.value);
            }
        }
        else {
            file << FormatValue(param.value);
        }
        if (param.comment) file << " # " << *param.comment;
        file << '\n';
    });
}
