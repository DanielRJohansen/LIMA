#include "Programs.h"

#include "Filehandling.h"

#include <glm.hpp>

#include <algorithm>
#include <array>
#include <cctype>
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

struct TerminalReplacement {
    std::string sourceName;
    std::string targetName;
    std::string type;
    double charge{};
    double mass{};
};

struct BondedTypeDefaults {
    int bond{};
    int angle{};
    int properDihedral{};
    int improperDihedral{};
    int generateAllDihedrals{};
    int exclusions{};
    int generateHydrogenPairs{};
    int removeDihedralsWithImpropers{};
};

BondedTypeDefaults readBondedTypeDefaults(const fs::path& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error(std::format("Failed to open residue database {}", path.string()));
    bool inBondedTypes = false;
    for (std::string line; std::getline(input, line);) {
        line = withoutComment(std::move(line));
        if (line.empty()) continue;
        if (line.front() == '[' && line.back() == ']') {
            inBondedTypes = trim(std::string_view(line).substr(1, line.size() - 2)) == "bondedtypes";
            continue;
        }
        if (inBondedTypes) {
            const auto fields = words(line);
            if (fields.size() < 8) break;
            return { std::stoi(fields[0]), std::stoi(fields[1]), std::stoi(fields[2]), std::stoi(fields[3]),
                std::stoi(fields[4]), std::stoi(fields[5]), std::stoi(fields[6]), std::stoi(fields[7]) };
        }
    }
    throw std::runtime_error(std::format("No [ bondedtypes ] defaults in {}", path.string()));
}

struct TerminalAddition {
    HydrogenInstruction instruction;
    std::string type;
    double charge{};
    int chargeGroup{};
    double mass{};
};

struct TerminalPatch {
    std::vector<TerminalReplacement> replacements;
    std::vector<TerminalAddition> additions;
    std::set<std::string> deletions;
    std::vector<NamedInteraction> bonds;
    std::vector<NamedInteraction> impropers;
    std::vector<NamedInteraction> cmaps;
};

std::string lowercase(std::string value) {
    std::ranges::transform(value, value.begin(), [](const unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return value;
}

std::unordered_map<std::string, TerminalPatch> readTerminalPatches(const fs::path& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error(std::format("Failed to open terminal database {}", path.string()));

    std::unordered_map<std::string, TerminalPatch> result;
    std::string patchName;
    std::string subsection;
    std::optional<HydrogenInstruction> pendingAddition;
    for (std::string line; std::getline(input, line);) {
        line = withoutComment(std::move(line));
        if (line.empty()) continue;
        if (line.front() == '[' && line.back() == ']') {
            const std::string section = trim(std::string_view(line).substr(1, line.size() - 2));
            const std::string normalized = lowercase(section);
            if (normalized == "replace" || normalized == "add" || normalized == "delete"
                || normalized == "bonds" || normalized == "impropers" || normalized == "cmap") {
                subsection = normalized;
            }
            else {
                patchName = section;
                subsection.clear();
                result.try_emplace(patchName);
            }
            pendingAddition.reset();
            continue;
        }
        if (patchName.empty() || subsection.empty()) continue;

        const auto fields = words(line);
        auto& patch = result.at(patchName);
        if (subsection == "replace" && fields.size() >= 4) {
            const bool renamed = fields.size() >= 5;
            const std::string& type = fields[renamed ? 2 : 1];
            const double declaredMass = std::stod(fields[renamed ? 3 : 2]);
            patch.replacements.push_back({ fields[0], renamed ? fields[1] : fields[0], type,
                std::stod(fields[renamed ? 4 : 3]), declaredMass });
        }
        else if (subsection == "add") {
            if (!pendingAddition) {
                if (fields.size() < 4) throw std::runtime_error(std::format("Malformed terminal add instruction: {}", line));
                HydrogenInstruction instruction;
                instruction.count = std::stoi(fields[0]);
                instruction.type = std::stoi(fields[1]);
                instruction.namePrefix = fields[2];
                instruction.controls.assign(fields.begin() + 3, fields.end());
                pendingAddition = std::move(instruction);
            }
            else {
                if (fields.size() < 4) throw std::runtime_error(std::format("Malformed terminal atom definition: {}", line));
                patch.additions.push_back({ *pendingAddition, fields[0], std::stod(fields[2]), std::stoi(fields[3]),
                    std::stod(fields[1]) });
                pendingAddition.reset();
            }
        }
        else if (subsection == "delete" && !fields.empty()) patch.deletions.insert(fields[0]);
        else if (subsection == "bonds" && fields.size() >= 2) patch.bonds.push_back({ fields[0], fields[1] });
        else if (subsection == "impropers" && fields.size() >= 4) patch.impropers.push_back({ fields[0], fields[1], fields[2], fields[3] });
        else if (subsection == "cmap" && fields.size() >= 5) patch.cmaps.push_back({ fields[0], fields[1], fields[2], fields[3], fields[4] });
    }
    return result;
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
    std::vector<HydrogenInstruction> terminalCoordinateInstructions;
    std::vector<NamedInteraction> terminalBonds;
    std::vector<NamedInteraction> terminalImpropers;
    std::vector<NamedInteraction> terminalCmaps;
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

void applyTerminalPatch(std::vector<TemplateAtom>& atoms, OutputResidue& output, const TerminalPatch& patch) {
    std::vector<TemplateAtom> modified;
    modified.reserve(atoms.size());
    for (auto atom : atoms) {
        if (patch.deletions.contains(atom.name)) continue;
        if (const auto replacement = std::ranges::find(patch.replacements, atom.name, &TerminalReplacement::sourceName);
            replacement != patch.replacements.end()) {
            atom.name = replacement->targetName;
            atom.type = replacement->type;
            atom.charge = replacement->charge;
            atom.mass = replacement->mass;
        }
        modified.push_back(std::move(atom));
    }

    for (const auto& addition : patch.additions) {
        const auto names = generatedNames(addition.instruction);
        std::size_t insertion = 0;
        if (!addition.instruction.controls.empty()) {
            const auto parent = std::ranges::find(modified, addition.instruction.controls.front(), &TemplateAtom::name);
            if (parent != modified.end()) insertion = static_cast<std::size_t>(std::distance(modified.begin(), parent) + 1);
        }
        for (const auto& name : names) {
            const auto existing = std::ranges::find(modified, name, &TemplateAtom::name);
            if (existing != modified.end()) insertion = std::max(insertion, static_cast<std::size_t>(std::distance(modified.begin(), existing) + 1));
        }
        for (const auto& name : names) {
            if (std::ranges::find(modified, name, &TemplateAtom::name) == modified.end()) {
                modified.insert(modified.begin() + static_cast<std::ptrdiff_t>(insertion++),
                    TemplateAtom{ name, addition.type, addition.charge, addition.chargeGroup, addition.mass });
            }
            if (!addition.instruction.controls.empty()) {
                output.terminalBonds.push_back({ addition.instruction.controls.front(), name });
            }
        }
        output.terminalCoordinateInstructions.push_back(addition.instruction);
    }
    output.terminalBonds.insert(output.terminalBonds.end(), patch.bonds.begin(), patch.bonds.end());
    output.terminalImpropers.insert(output.terminalImpropers.end(), patch.impropers.begin(), patch.impropers.end());
    output.terminalCmaps.insert(output.terminalCmaps.end(), patch.cmaps.begin(), patch.cmaps.end());
    atoms = std::move(modified);
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
            instructions.insert(instructions.begin(), residue.terminalCoordinateInstructions.begin(),
                residue.terminalCoordinateInstructions.end());

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
    const std::unordered_map<std::string, std::vector<HydrogenInstruction>>& hydrogenDatabase,
    const std::unordered_map<std::string, TerminalPatch>& nTerminalPatches,
    const std::unordered_map<std::string, TerminalPatch>& cTerminalPatches) {
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
        OutputResidue output;
        output.sourceName = source.name;
        output.rtpName = rtpName;
        if (residueIndex == 0) {
            const std::string patchName = source.name == "GLY" ? "GLY-NH3+" : source.name == "PRO" ? "PRO-NH2+" : "NH3+";
            applyTerminalPatch(definitions, output, nTerminalPatches.at(patchName));
        }
        if (residueIndex == static_cast<int>(pdb.residues.size()) - 1) {
            applyTerminalPatch(definitions, output, cTerminalPatches.at("COO-"));
        }

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
        for (const auto& names : residues[residue].terminalBonds) {
            if (auto bond = resolveInteraction<2>(names, residues, offsets, residue)) {
                if ((*bond)[1] < (*bond)[0]) std::swap((*bond)[0], (*bond)[1]);
                bondSet.insert(*bond);
            }
        }
        for (const auto& names : residues[residue].terminalImpropers) {
            if (auto improper = resolveInteraction<4>(names, residues, offsets, residue)) result.impropers.push_back(*improper);
        }
        for (const auto& names : residues[residue].terminalCmaps) {
            if (auto cmap = resolveInteraction<5>(names, residues, offsets, residue)) result.cmaps.push_back(*cmap);
        }
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

std::string_view waterTopologyFilename(const Programs::WaterModel model) {
    switch (model) {
        case Programs::WaterModel::Tip3p: return "tip3p.itp";
        case Programs::WaterModel::Tip4p: return "tip4p.itp";
        case Programs::WaterModel::Tips3p: return "tips3p.itp";
        case Programs::WaterModel::Tip5p: return "tip5p.itp";
        case Programs::WaterModel::Spc: return "spc.itp";
        case Programs::WaterModel::Spce: return "spce.itp";
    }
    throw std::runtime_error("Unknown water model");
}

void writeTopology(const fs::path& path, const fs::path& positionRestraints,
                   const std::string& title, char chain,
                   const std::vector<OutputResidue>& residues, const Interactions& interactions,
                   Programs::WaterModel waterModel, const BondedTypeDefaults& bondedTypes) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error(std::format("Failed to create {}", path.string()));
    const std::string moleculeName = chain == ' ' ? "Protein" : std::string("Protein_chain_") + chain;
    output << "; Topology generated by LIMA pdb2gmx\n\n"
              "; Include forcefield parameters\n"
              "#include \"charmm27.ff/forcefield.itp\"\n\n"
              "[ moleculetype ]\n"
              "; Name            nrexcl\n"
           << moleculeName << "     " << bondedTypes.exclusions << "\n\n"
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
    writeInteractionSection(output, "bonds", interactions.bonds, bondedTypes.bond);
    writeInteractionSection(output, "pairs", interactions.pairs, 1);
    writeInteractionSection(output, "angles", interactions.angles, bondedTypes.angle);
    writeInteractionSection(output, "dihedrals", interactions.propers, bondedTypes.properDihedral);
    writeInteractionSection(output, "dihedrals", interactions.impropers, bondedTypes.improperDihedral);
    writeInteractionSection(output, "cmap", interactions.cmaps, 1);
    output << "#ifdef POSRES\n#include \"" << positionRestraints.filename().string() << "\"\n#endif\n\n"
              "; Include water topology\n#include \"charmm27.ff/" << waterTopologyFilename(waterModel) << "\"\n\n"
              "; Include topology for ions\n#include \"charmm27.ff/ions.itp\"\n\n"
              "[ system ]\n; Name\n" << title << "\n\n"
              "[ molecules ]\n; Compound        #mols\n" << moleculeName << "     1\n";
}

} // namespace

void Programs::pdb2gmx(const fs::path& pdbfile, std::optional<std::string> name, WaterModel waterModel) {
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
    const auto bondedTypes = readBondedTypeDefaults(forcefieldDirectory / "aminoacids.rtp");
    const auto templates = readResidueTemplates(forcefieldDirectory / "aminoacids.rtp", masses);
    const auto hydrogenDatabase = readHydrogenDatabase(forcefieldDirectory / "aminoacids.hdb");
    const auto nTerminalPatches = readTerminalPatches(forcefieldDirectory / "aminoacids.n.tdb");
    const auto cTerminalPatches = readTerminalPatches(forcefieldDirectory / "aminoacids.c.tdb");
    const PdbInput pdb = readPdb(pdbfile);
    const auto outputResidues = makeAtoms(pdb, templates, hydrogenDatabase, nTerminalPatches, cTerminalPatches);
    const auto interactions = makeInteractions(outputResidues, templates);

    const fs::path directory = pdbfile.parent_path().empty() ? fs::current_path() : pdbfile.parent_path();
    const std::string basename = name.value_or("");
    const fs::path groPath = directory / (basename.empty() ? "conf.gro" : basename + ".gro");
    const fs::path topPath = directory / (basename.empty() ? "topol.top" : basename + ".top");
    const fs::path posrePath = directory / (basename.empty() ? "posre.itp" : basename + "_posre.itp");

    writeGro(groPath, pdb.title, outputResidues, pdb.box);
    writePositionRestraints(posrePath, outputResidues);
    writeTopology(topPath, posrePath, pdb.title, pdb.residues.front().chain, outputResidues, interactions,
        waterModel, bondedTypes);
}
