#include "SimParams.h"
#include "Format.h"

#include "Filehandling.h"

#include <algorithm>
#include <array>
#include <charconv>
#include <cctype>
#include <concepts>
#include <format>
#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>

namespace fs = std::filesystem;
using Dictionary = std::unordered_map<std::string, std::string>;

SimParams SimParams::BasicEMSimParams(float forceTolerance) {
    SimParams params{};
    params.em_variant = true;
    params.dt = 1.5f * FEMTO_TO_NANO;
    params.emtol = forceTolerance;
    params.nstenergy = 50;
    params.enable_electrostatics = true;
    params.nsteps = 20000;
    return params;
}

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
        throw std::runtime_error(Lima::Format("Invalid value for '{}': '{}'", key, input));
    return value;
}

template<typename T>
void ParseValue(const std::string_view key, const std::string_view input, T& value) {
    if constexpr (std::same_as<T, bool>) {
        if (input == "true") value = true;
        else if (input == "false") value = false;
        else throw std::runtime_error(Lima::Format("Invalid boolean for '{}': '{}'", key, input));
    }
    else if constexpr (std::integral<T> || std::floating_point<T>) {
        value = ParseNumber<T>(key, input);
    }
    else if constexpr (std::same_as<T, BoundaryConditionSelect>) {
        if (input == "pbc") value = PBC;
        else if (input == "nobc") value = NoBC;
        else throw std::runtime_error(Lima::Format("Invalid boundary condition: '{}'", input));
    }
    else if constexpr (std::same_as<T, ConstraintsSelect>) {
        if (input == "none") value = ConstraintsSelect::None;
        else if (input == "h-bonds") value = ConstraintsSelect::HBonds;
        else throw std::runtime_error(Lima::Format("Invalid constraints: '{}', expected none or h-bonds", input));
    }
    else if constexpr (std::same_as<T, TcouplSelect>) {
        if (input == "no" || input == "none") value = TcouplSelect::No;
        else if (input == "v-rescale") value = TcouplSelect::VRescale;
        else throw std::runtime_error(Lima::Format("Unsupported tcoupl: '{}', LIMA has no or v-rescale", input));
    }
    else if constexpr (std::same_as<T, PcouplSelect>) {
        if (input == "no" || input == "none") value = PcouplSelect::No;
        else if (input == "berendsen") value = PcouplSelect::Berendsen;
        else throw std::runtime_error(Lima::Format("Unsupported pcoupl: '{}', LIMA has no or berendsen (isotropic)", input));
    }
    else if constexpr (std::same_as<T, ColoringMethod>) {
        static const std::unordered_map<std::string_view, ColoringMethod> values{
            { "atomname", ColoringMethod::Atomname },
            { "charge", ColoringMethod::Charge },
            { "gradientfromatomid", ColoringMethod::GradientFromAtomid },
            { "persistentclusterid", ColoringMethod::PersistentClusterId },
            { "forcemagnitude", ColoringMethod::ForceMagnitude },
            { "newcartoon", ColoringMethod::NewCartoon },
            { "lipids", ColoringMethod::Lipids },
            { "licorice", ColoringMethod::Licorice },
            { "surface", ColoringMethod::Surface }
        };
        const auto found = values.find(input);
        if (found == values.end()) throw std::runtime_error(Lima::Format("Invalid coloring method: '{}'", input));
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
            else throw std::runtime_error(Lima::Format("Invalid supernatural force: '{}'", item));
        }
    }
}

std::string FormatValue(const bool value) { return value ? "true" : "false"; }
std::string FormatValue(const BoundaryConditionSelect value) { return value == PBC ? "PBC" : "NoBC"; }

std::string FormatValue(const ConstraintsSelect value) { return value == ConstraintsSelect::HBonds ? "h-bonds" : "none"; }
std::string FormatValue(const TcouplSelect value) { return value == TcouplSelect::VRescale ? "v-rescale" : "no"; }
std::string FormatValue(const PcouplSelect value) { return value == PcouplSelect::Berendsen ? "berendsen" : "no"; }

std::string FormatValue(const ColoringMethod value) {
    switch (value) {
    case ColoringMethod::Atomname: return "Atomname";
    case ColoringMethod::Charge: return "Charge";
    case ColoringMethod::GradientFromAtomid: return "GradientFromAtomid";
    case ColoringMethod::PersistentClusterId: return "PersistentClusterId";
    case ColoringMethod::ForceMagnitude: return "ForceMagnitude";
    case ColoringMethod::NewCartoon: return "NewCartoon";
    case ColoringMethod::Lipids: return "Lipids";
    case ColoringMethod::Licorice: return "Licorice";
    case ColoringMethod::Surface: return "Surface";
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
    case SimParamSection::Barostat: return "Barostat parameters";
    case SimParamSection::Output: return "Output parameters";
    case SimParamSection::Analysis: return "Analysis";
    case SimParamSection::Debug: return "Debug parameters";
    }
    throw std::runtime_error("Invalid simulation parameter section");
}

template<typename T>
std::string FormatValue(const T value) {
    return Lima::Format("{}", value);
}

template<typename Function>
void ForEachParam(SimParams& params, Function&& function) {
    std::apply([&](auto&... param) { (function(param), ...); }, params.Params());
}

template<typename Function>
void ForEachParam(const SimParams& params, Function&& function) {
    std::apply([&](const auto&... param) { (function(param), ...); }, params.Params());
}

// An .mdp file as GROMACS reads it: ';' starts a comment, keys are case-insensitive with '-' and '_' equal (normalized to '-'),
// and values keep their inner spaces, as some have one entry per group. Values are lowercased too, as LIMA compares them so
Dictionary ReadMdp(const fs::path& path) {
    std::ifstream file{ path };
    if (!file)
        throw std::runtime_error("Failed to open file " + path.string());
    const auto Trim = [](std::string_view sv) {
        const size_t begin = sv.find_first_not_of(" \t\r");
        if (begin == std::string_view::npos) return std::string_view{};
        return sv.substr(begin, sv.find_last_not_of(" \t\r") - begin + 1);
    };
    Dictionary mdp;
    std::string line;
    while (std::getline(file, line)) {
        std::string_view sv{ line };
        sv = sv.substr(0, sv.find(';'));
        const size_t equals = sv.find('=');
        if (equals == std::string_view::npos)
            continue;
        std::string key = Lowercase(Trim(sv.substr(0, equals)));
        std::ranges::replace(key, '_', '-');
        std::string value;
        for (const char c : Trim(sv.substr(equals + 1))) {
            const bool space = c == ' ' || c == '\t';
            if (!space || (!value.empty() && value.back() != ' '))
                value += space ? ' ' : static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        }
        if (key.empty())
            continue;
        if (!mdp.emplace(key, value).second)
            throw std::runtime_error(Lima::Format("'{}' is given more than once in {}", key, path.string()));
    }
    return mdp;
}

void ParseMdp(const Dictionary& mdp, SimParams& params) {
    const auto parseIfPresent = [&mdp]<typename T>(const std::string_view key, T& target) {
        if (const auto found = mdp.find(std::string{ key }); found != mdp.end())
            ParseValue(key, found->second, target);
    };

    bool velocityVerlet = false;
    if (const auto found = mdp.find("integrator"); found != mdp.end()) {
        if (found->second == "md-vv") {
            params.em_variant.value = false;
            velocityVerlet = true;
        }
        else if (found->second == "md")
            throw std::runtime_error("Unsupported integrator: md is GROMACS leap-frog; LIMA uses velocity Verlet, set integrator = md-vv");
        else if (found->second == "steep") params.em_variant.value = true;
        else throw std::runtime_error("Unsupported integrator: " + found->second);
    }

    if (const auto found = mdp.find("nsteps"); found != mdp.end()) {
        // GROMACS' -1 runs without end. LIMA sizes its output buffers by the step count, so it needs one
        if (found->second == "-1")
            std::cerr << Lima::Format("Warning: nsteps = -1 (no limit) is not supported, running the default {} steps. Set "
                "nsteps, eg. with --set nsteps=...\n", params.nsteps.value);
        else
            ParseValue("nsteps", found->second, params.nsteps.value);
    }
    if (const auto found = mdp.find("dt"); found != mdp.end())
        params.dt.value = ParseNumber<float>("dt", found->second) * PICO_TO_NANO;
    parseIfPresent("emtol", params.emtol.value);
    parseIfPresent("nstlist", params.nstlist.value);
    parseIfPresent("verlet-buffer-tolerance", params.verlet_buffer_tolerance.value);
    parseIfPresent("rlist", params.rlist.value);
    parseIfPresent("constraints", params.constraints.value);

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
    parseIfPresent("fourierspacing", params.fourier_spacing.value);
    // GROMACS' multiple time stepping, with the PME mesh as its only level-2 force
    if (const auto found = mdp.find("mts"); found != mdp.end()) {
        if (found->second == "yes") params.mts.value = true;
        else if (found->second == "no") params.mts.value = false;
        else throw std::runtime_error("Invalid value for 'mts': " + found->second);
    }
    if (params.mts) {
        if (const auto found = mdp.find("mts-levels"); found != mdp.end() && found->second != "2")
            throw std::runtime_error("Unsupported mts-levels: " + found->second + ", LIMA has 2");
        if (const auto found = mdp.find("mts-level2-forces"); found != mdp.end() && found->second != "longrange-nonbonded")
            throw std::runtime_error("Unsupported mts-level2-forces: " + found->second + ", LIMA has longrange-nonbonded");
        parseIfPresent("mts-level2-factor", params.mts_level2_factor.value);
    }
    // LIMA has one cutoff for LJ and Coulomb
    if (const auto found = mdp.find("rvdw"); found != mdp.end())
        if (ParseNumber<float>("rvdw", found->second) != params.cutoff_nm)
            throw std::runtime_error(Lima::Format("rvdw = {} differs from rcoulomb = {}. LIMA uses one cutoff for both", found->second, params.cutoff_nm.value));

    // The trajectory file has its own interval. 0 means none in GROMACS; LIMA writes one only for mdrun --trajectory
    for (const std::string_view key : { "nstxout", "nstxout-compressed" })
        if (const auto found = mdp.find(std::string{ key }); found != mdp.end())
            if (const int interval = ParseNumber<int>(key, found->second); interval > 0)
                params.nstxout_compressed.value = interval;	// nstxout-compressed wins, as it comes last
    constexpr std::array outputKeys{ "nstvout", "nstenergy", "nstlog", "nstcomm" };
    std::optional<int> outputInterval;
    for (const std::string_view key : outputKeys) {
        if (const auto found = mdp.find(std::string{ key }); found != mdp.end()) {
            const int interval = ParseNumber<int>(key, found->second);
            if (outputInterval && *outputInterval != interval)
                throw std::runtime_error(Lima::Format("Output interval '{}' differs from earlier intervals", key));
            outputInterval = interval;
        }
    }
    if (outputInterval) params.nstenergy.value = *outputInterval;

    parseIfPresent("tcoupl", params.tcoupl.value);
    parseIfPresent("nsttcouple", params.nsttcouple.value);
    // GROMACS gives these per coupling group. LIMA couples the whole system, so the groups must agree
    const auto parseSameForAllGroups = [&mdp](const std::string_view key, float& target) {
        const auto found = mdp.find(std::string{ key });
        if (found == mdp.end())
            return;
        std::optional<float> value;
        std::string_view rest = found->second;
        while (!rest.empty()) {
            const size_t start = rest.find_first_not_of(" \t");
            if (start == std::string_view::npos)
                break;
            rest = rest.substr(start);
            const size_t end = std::min(rest.find_first_of(" \t"), rest.size());
            const float groupValue = ParseNumber<float>(key, rest.substr(0, end));
            if (value && *value != groupValue)
                throw std::runtime_error(Lima::Format("'{}' differs between coupling groups, LIMA couples the whole system", key));
            value = groupValue;
            rest = rest.substr(end);
        }
        if (value)
            target = *value;
    };
    parseSameForAllGroups("tau-t", params.tau_t.value);
    parseSameForAllGroups("ref-t", params.ref_t.value);
    // GROMACS' default -1 picks a random seed, which would make runs irreproducible, so it keeps LIMA's fixed seed
    const auto parseSeed = [&mdp](const std::string_view key, int64_t& target) {
        if (const auto found = mdp.find(std::string{ key }); found != mdp.end()) {
            const int64_t seed = ParseNumber<int64_t>(key, found->second);
            if (seed != -1)
                target = seed;
        }
    };
    parseSeed("ld-seed", params.ld_seed.value);
    parseSeed("gen-seed", params.gen_seed.value);

    parseIfPresent("pcoupl", params.pcoupl.value);
    if (const auto found = mdp.find("pcoupltype"); found != mdp.end() && Lowercase(found->second) != "isotropic")
        throw std::runtime_error("Unsupported pcoupltype: " + found->second + ". LIMA has isotropic pressure coupling");
    parseIfPresent("nstpcouple", params.nstpcouple.value);
    parseSameForAllGroups("tau-p", params.tau_p.value);
    parseSameForAllGroups("ref-p", params.ref_p.value);
    parseSameForAllGroups("compressibility", params.compressibility.value);

    // GROMACS couples every step with velocity Verlet, irrespective of nsttcouple.
    if (velocityVerlet && params.tcoupl != TcouplSelect::No)
        params.nsttcouple.value = 1;
}

} // namespace

SimParams::SimParams(const fs::path& path) {
    if (path.extension() == ".mdp") {
        ParseMdp(ReadMdp(path), *this);
        return;
    }
    const Dictionary dictionary = FileUtils::parseINIFile(path.string(), true);

    for (const auto& [key, value] : dictionary)
        Set(key, value);
}

void SimParams::Set(std::string_view key, std::string_view value) {
    // '-' and '_' are equal in keys, as in GROMACS' .mdp files, eg. nstxout-compressed
    std::string normalizedKey = Lowercase(key);
    std::ranges::replace(normalizedKey, '-', '_');
    // Values are compared lowercase, as DumpToFile writes eg. PBC and Atomname, and --set passes them on as typed
    const std::string lowercaseValue = Lowercase(value);
    bool known = false;
    ForEachParam(*this, [&](auto& param) {
        if (known || Lowercase(param.name) != normalizedKey) return;
        known = true;
        if constexpr (std::same_as<typename std::remove_cvref_t<decltype(param)>::Type, float>) {
            if (param.name == "dt") {
                // [ps], as GROMACS. Via fs, so eg. 0.002 gives the same float as 2 fs did, when dt was in fs
                const double picoseconds = ParseNumber<double>(param.name, lowercaseValue);
                if (picoseconds > 0.1)
                    throw std::runtime_error(Lima::Format("dt = {} ps is too large. dt is in ps, as in GROMACS: eg. 0.002 for 2 fs", value));
                param.value = static_cast<float>(picoseconds * 1000.0) * FEMTO_TO_NANO;
                return;
            }
        }
        ParseValue(param.name, lowercaseValue, param.value);
    });
    if (!known)
        throw std::runtime_error(Lima::Format("Unknown simulation parameter '{}={}'", key, value));
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
                // Via fs, which rounds to eg. 2 exactly, so the file reads 0.002 and not 0.0019999999
                const float femtoseconds = static_cast<float>(param.value * NANO_TO_FEMTO);
                file << FormatValue(static_cast<float>(femtoseconds / 1000.0));
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
