#include "Programs.h"

#include "Filehandling.h"

#include <glm.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <format>
#include <iomanip>
#include <numbers>
#include <optional>
#include <ranges>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace {
std::string trim(std::string_view text) {
    const auto first = text.find_first_not_of(" \t\r\n");
    if (first == std::string_view::npos) return {};
    const auto last = text.find_last_not_of(" \t\r\n");
    return std::string(text.substr(first, last - first + 1));
}

std::vector<std::string> words(std::string_view line) {
    std::istringstream input{ std::string(line) };
    std::vector<std::string> result;
    for (std::string word; input >> word;) result.push_back(std::move(word));
    return result;
}

std::string withoutComment(std::string line) {
    if (const auto comment = line.find(';'); comment != std::string::npos) line.erase(comment);
    return trim(line);
}

struct PdbAtom {
    std::string name;
    std::string residueName;
    char chain{};
    int residueNumber{};
    char insertionCode{};
    glm::vec3 position;
};

struct PdbResidue {
    std::string name;
    char chain{};
    int number{};
    char insertionCode{};
    std::unordered_map<std::string, PdbAtom> atoms;
};

struct CrystalBox { std::array<double, 9> groValues{}; };

struct PdbInput {
    std::string title;
    std::vector<PdbResidue> residues;
    std::optional<CrystalBox> box;
};

PdbInput readPdb(const fs::path& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error(std::format("Failed to open PDB file {}", path.string()));

    PdbInput result;
    for (std::string line; std::getline(input, line);) {
        if (line.starts_with("COMPND") && line.size() > 20) {
            const std::string contents = trim(std::string_view(line).substr(10));
            if (const auto marker = contents.find("MOLECULE:"); marker != std::string::npos && result.title.empty()) {
                result.title = trim(contents.substr(marker + 9));
                if (!result.title.empty() && result.title.back() == ';') result.title.pop_back();
            }
        }
        if (line.starts_with("CRYST1") && line.size() >= 54) {
            const double a = std::stod(line.substr(6, 9)) * 0.1;
            const double b = std::stod(line.substr(15, 9)) * 0.1;
            const double c = std::stod(line.substr(24, 9)) * 0.1;
            const double alpha = std::stod(line.substr(33, 7)) * std::numbers::pi / 180.0;
            const double beta = std::stod(line.substr(40, 7)) * std::numbers::pi / 180.0;
            const double gamma = std::stod(line.substr(47, 7)) * std::numbers::pi / 180.0;
            const double bx = b * std::cos(gamma);
            const double by = b * std::sin(gamma);
            const double cx = c * std::cos(beta);
            const double cy = c * (std::cos(alpha) - std::cos(beta) * std::cos(gamma)) / std::sin(gamma);
            const double cz = std::sqrt(std::max(0.0, c * c - cx * cx - cy * cy));
            result.box = CrystalBox{ { a, by, cz, 0.0, 0.0, bx, 0.0, cx, cy } };
        }

        if (!line.starts_with("ATOM  ") || line.size() < 54) continue;
        const char alternateLocation = line.size() > 16 ? line[16] : ' ';
        if (alternateLocation != ' ' && alternateLocation != 'A') continue;

        PdbAtom atom;
        atom.name = trim(std::string_view(line).substr(12, 4));
        if (atom.name == "H") atom.name = "HN"; // CHARMM27 aminoacids.arn mapping.
        atom.residueName = trim(std::string_view(line).substr(17, 3));
        // Standard PDB and CHARMM use different names for isoleucine's terminal carbon.
        if (atom.residueName == "ILE" && atom.name == "CD1") atom.name = "CD";
        atom.chain = line.size() > 21 ? line[21] : ' ';
        atom.residueNumber = std::stoi(line.substr(22, 4));
        atom.insertionCode = line.size() > 26 ? line[26] : ' ';
        atom.position = { std::stof(line.substr(30, 8)) * 0.1f,
                          std::stof(line.substr(38, 8)) * 0.1f,
                          std::stof(line.substr(46, 8)) * 0.1f };

        if (result.residues.empty()
            || result.residues.back().chain != atom.chain
            || result.residues.back().number != atom.residueNumber
            || result.residues.back().insertionCode != atom.insertionCode) {
            result.residues.push_back(PdbResidue{ atom.residueName, atom.chain, atom.residueNumber, atom.insertionCode, {} });
        }
        result.residues.back().atoms.try_emplace(atom.name, std::move(atom));
    }

    if (result.residues.empty()) throw std::runtime_error(std::format("No ATOM records were found in {}", path.string()));
    if (result.title.empty()) result.title = path.stem().string();
    return result;
}

struct TemplateAtom {
    std::string name;
    std::string type;
    double charge{};
    int chargeGroup{};
    double mass{};
};

using NamedInteraction = std::vector<std::string>;

struct ResidueTemplate {
    std::vector<TemplateAtom> atoms;
    std::vector<NamedInteraction> bonds;
    std::vector<NamedInteraction> impropers;
    std::vector<NamedInteraction> cmaps;
};

std::unordered_map<std::string, double> readAtomMasses(const fs::path& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error(std::format("Failed to open atom type database {}", path.string()));
    std::unordered_map<std::string, double> masses;
    for (std::string line; std::getline(input, line);) {
        const auto fields = words(withoutComment(std::move(line)));
        if (fields.size() >= 2) {
            try { masses[fields[0]] = std::stod(fields[1]); }
            catch (const std::invalid_argument&) { /* Header line. */ }
        }
    }
    return masses;
}

std::unordered_map<std::string, ResidueTemplate> readResidueTemplates(
    const fs::path& path, const std::unordered_map<std::string, double>& masses) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error(std::format("Failed to open residue database {}", path.string()));

    std::unordered_map<std::string, ResidueTemplate> result;
    std::string residue;
    std::string subsection;
    for (std::string line; std::getline(input, line);) {
        line = withoutComment(std::move(line));
        if (line.empty()) continue;
        if (line.front() == '[' && line.back() == ']') {
            const std::string section = trim(std::string_view(line).substr(1, line.size() - 2));
            if (section == "atoms" || section == "bonds" || section == "impropers" || section == "cmap"
                || section == "angles" || section == "dihedrals" || section == "exclusions") {
                subsection = section;
            }
            else if (section == "bondedtypes") {
                residue.clear();
                subsection.clear();
            }
            else {
                residue = section;
                subsection.clear();
                result.try_emplace(residue);
            }
            continue;
        }
        if (residue.empty()) continue;

        const auto fields = words(line);
        auto& target = result.at(residue);
        if (subsection == "atoms" && fields.size() >= 4) {
            const auto mass = masses.find(fields[1]);
            if (mass == masses.end()) throw std::runtime_error(std::format("No mass for CHARMM27 atom type {}", fields[1]));
            target.atoms.push_back({ fields[0], fields[1], std::stod(fields[2]), std::stoi(fields[3]), mass->second });
        }
        else if (subsection == "bonds" && fields.size() >= 2) target.bonds.push_back({ fields[0], fields[1] });
        else if (subsection == "impropers" && fields.size() >= 4) target.impropers.push_back({ fields[0], fields[1], fields[2], fields[3] });
        else if (subsection == "cmap" && fields.size() >= 5) target.cmaps.push_back({ fields[0], fields[1], fields[2], fields[3], fields[4] });
    }
    return result;
}

struct HydrogenInstruction {
    int count{};
    int type{};
    std::string namePrefix;
    std::vector<std::string> controls;
};

std::unordered_map<std::string, std::vector<HydrogenInstruction>> readHydrogenDatabase(const fs::path& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error(std::format("Failed to open hydrogen database {}", path.string()));
    std::vector<std::string> lines;
    for (std::string line; std::getline(input, line);) {
        line = withoutComment(std::move(line));
        if (!line.empty()) lines.push_back(std::move(line));
    }

    std::unordered_map<std::string, std::vector<HydrogenInstruction>> result;
    for (std::size_t i = 0; i < lines.size();) {
        const auto header = words(lines[i++]);
        if (header.size() < 2) continue;
        const int count = std::stoi(header[1]);
        auto& instructions = result[header[0]];
        for (int row = 0; row < count && i < lines.size(); ++row, ++i) {
            const auto fields = words(lines[i]);
            if (fields.size() < 4) throw std::runtime_error(std::format("Malformed hydrogen database line: {}", lines[i]));
            HydrogenInstruction instruction;
            instruction.count = std::stoi(fields[0]);
            instruction.type = std::stoi(fields[1]);
            instruction.namePrefix = fields[2];
            instruction.controls.assign(fields.begin() + 3, fields.end());
            instructions.push_back(std::move(instruction));
        }
    }
    return result;
}

std::vector<std::string> generatedNames(const HydrogenInstruction& instruction) {
    if (instruction.count == 1) return { instruction.namePrefix };
    std::vector<std::string> names;
    for (int i = 1; i <= instruction.count; ++i) names.push_back(instruction.namePrefix + std::to_string(i));
    return names;
}

// Geometry rules are the standard GROMACS hydrogen-database construction
// rules. Distances are in nm and match the CHARMM27 pdb2gmx defaults.
std::vector<glm::vec3> calculateAddedPositions(int type, const std::vector<glm::vec3>& controls, int outputCount) {
    constexpr float distanceH = 0.1f;
    const float tetrahedral = std::acos(-1.f / 3.f);
    const float planar = 2.f * std::numbers::pi_v<float> / 3.f;
    std::vector<glm::vec3> output(static_cast<std::size_t>(outputCount));

    if (type == 1) {
        
        const glm::vec3 direction = glm::normalize(glm::normalize(controls[0] - controls[1]) + glm::normalize(controls[0] - controls[2]));
        output[0] = controls[0] + direction * distanceH;
        return output;
    }
    if (type == 5) {
        const glm::vec3 center = (controls[1] + controls[2] + controls[3]) / 3.f;
        output[0] = controls[0] + glm::normalize(controls[0] - center) * distanceH;
        return output;
    }
    if (type == 6) {
        const glm::vec3 bisector = controls[0] - (controls[1] + controls[2]) * 0.5f;
        const glm::vec3 normal = glm::cross(controls[0] - controls[1], controls[0] - controls[2]);
        const glm::vec3 a = glm::normalize(bisector) * (distanceH * std::cos(tetrahedral / 2.f));
        const glm::vec3 b = glm::normalize(normal) * (distanceH * std::sin(tetrahedral / 2.f));
        output[0] = controls[0] + a + b;
        output[1] = controls[0] + a - b;
        return output;
    }

    if (type == 2 || type == 3 || type == 4 || type == 8) {
        const glm::vec3 axis = glm::normalize(controls[0] - controls[1]);
        const glm::vec3 planeNormal = glm::normalize(glm::cross(axis, controls[1] - controls[2]));
        const glm::vec3 perpendicular = glm::cross(planeNormal, axis);
        if (type == 2) {
            output[0] = controls[0] + perpendicular * (distanceH * std::sin(tetrahedral))
                      - axis * (distanceH * std::cos(tetrahedral));
        }
        else if (type == 3) {
            output[0] = controls[0] - perpendicular * (distanceH * std::sin(planar))
                      - axis * (distanceH * std::cos(planar));
            output[1] = controls[0] + perpendicular * (distanceH * std::sin(planar))
                      - axis * (distanceH * std::cos(planar));
        }
        else if (type == 4) {
            const float side = distanceH * std::sin(tetrahedral);
            const glm::vec3 along = axis * (-distanceH * std::cos(tetrahedral));
            output[0] = controls[0] + perpendicular * side + along;
            output[1] = controls[0] - perpendicular * (side * 0.5f) + planeNormal * (side * std::sqrt(3.0f) * 0.5f) + along;
            if (outputCount >= 3) {
                output[2] = controls[0] - perpendicular * (side * 0.5f) - planeNormal * (side * std::sqrt(3.0f) * 0.5f) + along;
            }
        }
        else {
            constexpr float distanceO = 0.136f;
            const float angle = 117.f * std::numbers::pi_v<float> / 180.f;
            output[0] = controls[0] - perpendicular * (distanceO * std::sin(angle))
                      - axis * (distanceO * std::cos(angle));
            output[1] = controls[0] + perpendicular * (distanceO * std::sin(angle))
                      - axis * (distanceO * std::cos(angle));
        }
        return output;
    }
    throw std::runtime_error(std::format("Unsupported CHARMM hydrogen construction type {}", type));
}

struct OutputAtom {
    int residue{};
    std::string residueName;
    std::string name;
    std::string type;
    double charge{};
    double mass{};
    std::optional<glm::vec3> position;
};

struct OutputResidue {
    std::string sourceName;
    std::string rtpName;
    std::vector<OutputAtom> atoms;
    std::unordered_map<std::string, int> localIndex;
};

std::string rtpNameFor(const std::string& pdbName) {
    if (pdbName == "HIS" || pdbName == "HISD" || pdbName == "HIS1") return "HSD";
    if (pdbName == "HISE") return "HSE";
    if (pdbName == "HISH") return "HSP";
    if (pdbName == "LYSN") return "LSN";
    if (pdbName == "ASPH") return "ASPP";
    if (pdbName == "GLUH") return "GLUP";
    return pdbName;
}

void applyNTerminus(std::vector<TemplateAtom>& atoms, const std::string& residueName) {
    const bool glycine = residueName == "GLY";
    const bool proline = residueName == "PRO";
    std::vector<TemplateAtom> modified;
    for (auto atom : atoms) {
        if (atom.name == "HN") continue;
        if (atom.name == "N") {
            if (proline) {
                atom.type = "NP";
                atom.charge = -0.07;
            }
            else {
                atom.type = "NH3";
                atom.charge = -0.3;
            }
            modified.push_back(atom);
            const int hydrogenCount = proline ? 2 : 3;
            const double hydrogenCharge = proline ? 0.24 : 0.33;
            for (int i = 1; i <= hydrogenCount; ++i) {
                modified.push_back({ "H" + std::to_string(i), "HC", hydrogenCharge, -1, 1.008 });
            }
            continue;
        }
        if (atom.name == "CA") {
            atom.type = glycine ? "CT2" : proline ? "CP1" : "CT1";
            atom.charge = glycine ? 0.13 : proline ? 0.16 : 0.21;
        }
        else if (atom.name == "HA" && !proline) {
            atom.type = "HB";
            atom.charge = 0.10;
        }
        else if (atom.name == "CD" && proline) {
            atom.type = "CP3";
            atom.charge = 0.16;
        }
        modified.push_back(std::move(atom));
    }
    atoms = std::move(modified);
}

void applyCTerminus(std::vector<TemplateAtom>& atoms) {
    for (auto& atom : atoms) {
        if (atom.name == "C") {
            atom.type = "CC";
            atom.charge = 0.34;
        }
        else if (atom.name == "O") {
            atom.name = "OT1";
            atom.type = "OC";
            atom.charge = -0.67;
            atom.mass = 15.9994;
        }
    }
    atoms.push_back({ "OT2", "OC", -0.67, -1, 15.9994 });
}

std::optional<glm::vec3> lookupPosition(const std::vector<OutputResidue>& residues, int residue, std::string name) {
    if (!name.empty() && (name.front() == '-' || name.front() == '+')) {
        residue += name.front() == '-' ? -1 : 1;
        name.erase(name.begin());
    }
    if (residue < 0 || residue >= static_cast<int>(residues.size())) return std::nullopt;
    if (name == "O" && residue == static_cast<int>(residues.size()) - 1) name = "OT1";
    const auto found = residues[residue].localIndex.find(name);
    if (found == residues[residue].localIndex.end()) return std::nullopt;
    return residues[residue].atoms[found->second].position;
}

void assignAddedCoordinates(
    std::vector<OutputResidue>& residues,
    const std::unordered_map<std::string, std::vector<HydrogenInstruction>>& hydrogenDatabase) {
    for (int pass = 0; pass < 4; ++pass) {
        bool changed = false;
        for (int residueIndex = 0; residueIndex < static_cast<int>(residues.size()); ++residueIndex) {
            auto& residue = residues[residueIndex];
            const auto databaseEntry = hydrogenDatabase.find(residue.rtpName);
            if (databaseEntry == hydrogenDatabase.end()) {
                throw std::runtime_error(std::format("No hydrogen database entry for residue {}", residue.rtpName));
            }
            std::vector<HydrogenInstruction> instructions = databaseEntry->second;
            if (residueIndex == 0) {
                instructions.insert(instructions.begin(), HydrogenInstruction{
                    residue.sourceName == "PRO" ? 2 : 3, 4, "H", { "N", "CA", "C" } });
            }

            for (const auto& instruction : instructions) {
                const auto names = generatedNames(instruction);
                bool needed = false;
                for (const auto& name : names) {
                    const auto atom = residue.localIndex.find(name);
                    needed = needed || (atom != residue.localIndex.end() && !residue.atoms[atom->second].position.has_value());
                }
                if (!needed) continue;

                std::vector<glm::vec3> controls;
                for (const auto& controlName : instruction.controls) {
                    const auto value = lookupPosition(residues, residueIndex, controlName);
                    if (!value) {
                        controls.clear();
                        break;
                    }
                    controls.push_back(*value);
                }
                if (controls.size() != instruction.controls.size()) continue;

                const auto positions = calculateAddedPositions(instruction.type, controls, instruction.count);
                for (int i = 0; i < instruction.count; ++i) {
                    const auto atom = residue.localIndex.find(names[i]);
                    if (atom != residue.localIndex.end() && !residue.atoms[atom->second].position) {
                        residue.atoms[atom->second].position = positions[i];
                        changed = true;
                    }
                }
            }
        }
        if (!changed) break;
    }

    auto& last = residues.back();
    if (const auto ot2 = last.localIndex.find("OT2"); ot2 != last.localIndex.end() && !last.atoms[ot2->second].position) {
        std::vector<glm::vec3> controls;
        for (const std::string name : { "C", "CA", "N" }) {
            const auto position = lookupPosition(residues, static_cast<int>(residues.size()) - 1, name);
            if (!position) throw std::runtime_error("Cannot construct the C-terminal OT2 atom");
            controls.push_back(*position);
        }
        last.atoms[ot2->second].position = calculateAddedPositions(8, controls, 2)[1];
    }

    for (const auto& residue : residues) {
        for (const auto& atom : residue.atoms) {
            if (!atom.position) {
                throw std::runtime_error(std::format(
                    "Atom {} is missing in residue {} {} and could not be constructed",
                    atom.name, residue.sourceName, atom.residue));
            }
        }
    }
}

std::vector<OutputResidue> makeAtoms(
    const PdbInput& pdb,
    const std::unordered_map<std::string, ResidueTemplate>& templates,
    const std::unordered_map<std::string, std::vector<HydrogenInstruction>>& hydrogenDatabase) {
    const char chain = pdb.residues.front().chain;
    if (std::ranges::any_of(pdb.residues, [chain](const PdbResidue& residue) { return residue.chain != chain; })) {
        throw std::runtime_error("pdb2gmx currently requires a single protein chain");
    }

    std::vector<OutputResidue> result;
    result.reserve(pdb.residues.size());
    for (int residueIndex = 0; residueIndex < static_cast<int>(pdb.residues.size()); ++residueIndex) {
        const auto& source = pdb.residues[residueIndex];
        const std::string rtpName = rtpNameFor(source.name);
        const auto templateEntry = templates.find(rtpName);
        if (templateEntry == templates.end()) {
            throw std::runtime_error(std::format("Residue {} has no CHARMM27 amino-acid template", source.name));
        }
        std::vector<TemplateAtom> definitions = templateEntry->second.atoms;
        if (residueIndex == 0) applyNTerminus(definitions, source.name);
        if (residueIndex == static_cast<int>(pdb.residues.size()) - 1) applyCTerminus(definitions);

        OutputResidue output{ source.name, rtpName, {}, {} };
        for (const auto& definition : definitions) {
            OutputAtom atom;
            atom.residue = residueIndex + 1;
            atom.residueName = source.name;
            atom.name = definition.name;
            atom.type = definition.type;
            atom.charge = definition.charge;
            atom.mass = definition.mass;
            const std::string sourceName = atom.name == "OT1" ? "O" : atom.name;
            if (const auto found = source.atoms.find(sourceName); found != source.atoms.end()) {
                atom.position = found->second.position;
            }
            else if (!atom.name.starts_with('H') && atom.name != "OT2") {
                throw std::runtime_error(std::format(
                    "Heavy atom {} is missing in residue {} {}", atom.name, source.name, source.number));
            }
            output.localIndex.emplace(atom.name, static_cast<int>(output.atoms.size()));
            output.atoms.push_back(std::move(atom));
        }
        result.push_back(std::move(output));
    }
    assignAddedCoordinates(result, hydrogenDatabase);
    return result;
}

using Bond = std::array<int, 2>;
using Angle = std::array<int, 3>;
using Dihedral = std::array<int, 4>;
using Cmap = std::array<int, 5>;

std::optional<int> resolveAtom(
    const std::vector<OutputResidue>& residues,
    const std::vector<int>& offsets,
    int residue,
    std::string name) {
    if (!name.empty() && (name.front() == '-' || name.front() == '+')) {
        residue += name.front() == '-' ? -1 : 1;
        name.erase(name.begin());
    }
    if (residue < 0 || residue >= static_cast<int>(residues.size())) return std::nullopt;
    if (name == "O" && residue == static_cast<int>(residues.size()) - 1) name = "OT1";
    const auto found = residues[residue].localIndex.find(name);
    if (found == residues[residue].localIndex.end()) return std::nullopt;
    return offsets[residue] + found->second + 1;
}

template<std::size_t N>
std::optional<std::array<int, N>> resolveInteraction(
    const NamedInteraction& names,
    const std::vector<OutputResidue>& residues,
    const std::vector<int>& offsets,
    int residue) {
    std::array<int, N> output{};
    for (std::size_t i = 0; i < N; ++i) {
        const auto atom = resolveAtom(residues, offsets, residue, names[i]);
        if (!atom) return std::nullopt;
        output[i] = *atom;
    }
    return output;
}

struct Interactions {
    std::vector<Bond> bonds;
    std::vector<Bond> pairs;
    std::vector<Angle> angles;
    std::vector<Dihedral> propers;
    std::vector<Dihedral> impropers;
    std::vector<Cmap> cmaps;
};

Interactions makeInteractions(
    const std::vector<OutputResidue>& residues,
    const std::unordered_map<std::string, ResidueTemplate>& templates) {
    std::vector<int> offsets(residues.size());
    int atomCount = 0;
    for (int i = 0; i < static_cast<int>(residues.size()); ++i) {
        offsets[i] = atomCount;
        atomCount += static_cast<int>(residues[i].atoms.size());
    }

    std::set<Bond> bondSet;
    Interactions result;
    for (int residue = 0; residue < static_cast<int>(residues.size()); ++residue) {
        const auto& residueTemplate = templates.at(residues[residue].rtpName);
        for (const auto& names : residueTemplate.bonds) {
            if (auto bond = resolveInteraction<2>(names, residues, offsets, residue)) {
                if ((*bond)[1] < (*bond)[0]) std::swap((*bond)[0], (*bond)[1]);
                bondSet.insert(*bond);
            }
        }
        for (const auto& names : residueTemplate.impropers) {
            if (auto improper = resolveInteraction<4>(names, residues, offsets, residue)) result.impropers.push_back(*improper);
        }
        for (const auto& names : residueTemplate.cmaps) {
            if (auto cmap = resolveInteraction<5>(names, residues, offsets, residue)) result.cmaps.push_back(*cmap);
        }
    }

    const auto firstN = resolveAtom(residues, offsets, 0, "N");
    for (int i = 1; i <= (residues.front().sourceName == "PRO" ? 2 : 3); ++i) {
        const auto hydrogen = resolveAtom(residues, offsets, 0, "H" + std::to_string(i));
        if (firstN && hydrogen) bondSet.insert({ *firstN, *hydrogen });
    }
    const int lastResidue = static_cast<int>(residues.size()) - 1;
    const auto carbon = resolveAtom(residues, offsets, lastResidue, "C");
    const auto ot2 = resolveAtom(residues, offsets, lastResidue, "OT2");
    if (carbon && ot2) {
        Bond bond{ *carbon, *ot2 };
        if (bond[1] < bond[0]) std::swap(bond[0], bond[1]);
        bondSet.insert(bond);
    }
    if (const auto ca = resolveAtom(residues, offsets, lastResidue, "CA"); carbon && ca && ot2) {
        const auto ot1 = resolveAtom(residues, offsets, lastResidue, "OT1");
        if (ot1) result.impropers.push_back({ *carbon, *ca, *ot2, *ot1 });
    }
    result.bonds.assign(bondSet.begin(), bondSet.end());

    std::vector<std::set<int>> neighbours(static_cast<std::size_t>(atomCount + 1));
    for (const auto& bond : result.bonds) {
        neighbours[bond[0]].insert(bond[1]);
        neighbours[bond[1]].insert(bond[0]);
    }

    std::set<Bond> oneThree;
    for (int center = 1; center <= atomCount; ++center) {
        const std::vector<int> adjacent(neighbours[center].begin(), neighbours[center].end());
        for (std::size_t i = 0; i < adjacent.size(); ++i) {
            for (std::size_t j = i + 1; j < adjacent.size(); ++j) {
                result.angles.push_back({ adjacent[i], center, adjacent[j] });
                oneThree.insert({ adjacent[i], adjacent[j] });
            }
        }
    }

    std::set<Dihedral> properSet;
    std::set<Bond> pairSet;
    for (const auto& central : result.bonds) {
        for (const int first : neighbours[central[0]]) {
            if (first == central[1]) continue;
            for (const int fourth : neighbours[central[1]]) {
                if (fourth == central[0] || fourth == first) continue;
                Dihedral value{ first, central[0], central[1], fourth };
                const Dihedral reverse{ fourth, central[1], central[0], first };
                if (reverse < value) value = reverse;
                properSet.insert(value);
                const Bond pair{ std::min(first, fourth), std::max(first, fourth) };
                if (!bondSet.contains(pair) && !oneThree.contains(pair)) pairSet.insert(pair);
            }
        }
    }
    result.propers.assign(properSet.begin(), properSet.end());
    result.pairs.assign(pairSet.begin(), pairSet.end());
    return result;
}

void writeGro(const fs::path& path, const std::string& title,
              const std::vector<OutputResidue>& residues, const std::optional<CrystalBox>& box) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error(std::format("Failed to create {}", path.string()));
    std::size_t atomCount = 0;
    for (const auto& residue : residues) atomCount += residue.atoms.size();
    output << title << '\n' << atomCount << '\n';
    int atomId = 1;
    for (const auto& residue : residues) {
        for (const auto& atom : residue.atoms) {
            output << std::setw(5) << std::right << (atom.residue % 100000)
                   << std::setw(5) << std::left << atom.residueName.substr(0, 5)
                   << std::setw(5) << std::right << atom.name.substr(0, 5)
                   << std::setw(5) << std::right << (atomId++ % 100000)
                   << std::setw(8) << std::fixed << std::setprecision(3) << atom.position->x
                   << std::setw(8) << atom.position->y
                   << std::setw(8) << atom.position->z << '\n';
        }
    }
    const auto values = box ? box->groValues : std::array<double, 9>{};
    const int valueCount = box ? 9 : 3;
    for (int i = 0; i < valueCount; ++i) output << std::setw(10) << std::fixed << std::setprecision(5) << values[i];
    output << '\n';
}

template<std::size_t N>
void writeInteractionSection(std::ofstream& output, std::string_view name,
                             const std::vector<std::array<int, N>>& interactions, int function) {
    output << "[ " << name << " ]\n";
    for (const auto& interaction : interactions) {
        for (const int atom : interaction) output << std::setw(6) << atom;
        output << std::setw(6) << function << " \n";
    }
    output << '\n';
}

void writePositionRestraints(const fs::path& path, const std::vector<OutputResidue>& residues) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error(std::format("Failed to create {}", path.string()));
    output << "; Position restraints for all non-hydrogen protein atoms.\n\n"
              "[ position_restraints ]\n"
              "; atom  type      fx      fy      fz\n";
    int atomId = 1;
    for (const auto& residue : residues) {
        for (const auto& atom : residue.atoms) {
            if (!atom.name.starts_with('H')) output << std::setw(6) << atomId << "     1  1000  1000  1000\n";
            ++atomId;
        }
    }
}

void writeTopology(const fs::path& path, const fs::path& positionRestraints,
                   const std::string& title, char chain,
                   const std::vector<OutputResidue>& residues, const Interactions& interactions) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error(std::format("Failed to create {}", path.string()));
    const std::string moleculeName = chain == ' ' ? "Protein" : std::string("Protein_chain_") + chain;
    output << "; Topology generated by LIMA pdb2gmx\n\n"
              "; Include forcefield parameters\n"
              "#include \"charmm27.ff/forcefield.itp\"\n\n"
              "[ moleculetype ]\n"
              "; Name            nrexcl\n"
           << moleculeName << "     3\n\n"
              "[ atoms ]\n"
              "; nr type resnr residue atom cgnr charge mass\n";
    int atomId = 1;
    double totalCharge = 0.0;
    for (const auto& residue : residues) {
        for (std::size_t i = 0; i < residue.atoms.size(); ++i) {
            const auto& atom = residue.atoms[i];
            totalCharge += atom.charge;
            output << std::setw(6) << atomId
                   << std::setw(11) << atom.type
                   << std::setw(7) << atom.residue
                   << std::setw(7) << atom.residueName
                   << std::setw(7) << atom.name
                   << std::setw(7) << atomId
                   << std::setw(11) << std::setprecision(6) << std::defaultfloat << atom.charge
                   << std::setw(11) << atom.mass;
            if (i + 1 == residue.atoms.size()) {
                output << "   ; qtot " << std::setprecision(0) << std::fixed << totalCharge << std::defaultfloat;
            }
            output << '\n';
            ++atomId;
        }
    }
    output << '\n';
    writeInteractionSection(output, "bonds", interactions.bonds, 1);
    writeInteractionSection(output, "pairs", interactions.pairs, 1);
    writeInteractionSection(output, "angles", interactions.angles, 5);
    writeInteractionSection(output, "dihedrals", interactions.propers, 9);
    writeInteractionSection(output, "dihedrals", interactions.impropers, 2);
    writeInteractionSection(output, "cmap", interactions.cmaps, 1);
    output << "#ifdef POSRES\n#include \"" << positionRestraints.filename().string() << "\"\n#endif\n\n"
              "; Include water topology\n#include \"charmm27.ff/tip3p.itp\"\n\n"
              "; Include topology for ions\n#include \"charmm27.ff/ions.itp\"\n\n"
              "[ system ]\n; Name\n" << title << "\n\n"
              "[ molecules ]\n; Compound        #mols\n" << moleculeName << "     1\n";
}

} // namespace

void Programs::pdb2gmx(const fs::path& pdbfile, std::optional<std::string> name) {
    if (pdbfile.extension() != ".pdb") {
        throw std::runtime_error(std::format("Expected a .pdb input file, got {}", pdbfile.string()));
    }
    if (!fs::is_regular_file(pdbfile)) {
        throw std::runtime_error(std::format("PDB input file does not exist: {}", pdbfile.string()));
    }
    if (name && (name->empty() || fs::path(*name).has_parent_path())) {
        throw std::runtime_error("pdb2gmx output name must be a non-empty basename");
    }

    const fs::path forcefieldDirectory = FileUtils::GetLimaDir() / "resources" / "forcefields" / "charmm27.ff";
    const auto masses = readAtomMasses(forcefieldDirectory / "atomtypes.atp");
    const auto templates = readResidueTemplates(forcefieldDirectory / "aminoacids.rtp", masses);
    const auto hydrogenDatabase = readHydrogenDatabase(forcefieldDirectory / "aminoacids.hdb");
    const PdbInput pdb = readPdb(pdbfile);
    const auto outputResidues = makeAtoms(pdb, templates, hydrogenDatabase);
    const auto interactions = makeInteractions(outputResidues, templates);

    const fs::path directory = pdbfile.parent_path().empty() ? fs::current_path() : pdbfile.parent_path();
    const std::string basename = name.value_or("");
    const fs::path groPath = directory / (basename.empty() ? "conf.gro" : basename + ".gro");
    const fs::path topPath = directory / (basename.empty() ? "topol.top" : basename + ".top");
    const fs::path posrePath = directory / (basename.empty() ? "posre.itp" : basename + "_posre.itp");

    writeGro(groPath, pdb.title, outputResidues, pdb.box);
    writePositionRestraints(posrePath, outputResidues);
    writeTopology(topPath, posrePath, pdb.title, pdb.residues.front().chain, outputResidues, interactions);
}
