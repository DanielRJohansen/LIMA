#pragma once

#include "Programs.h"
#include "argparser.h"

#include <format>
#include <stdexcept>

inline int togmx(int argc, char** argv) {
    const std::string helpText = R"(
Usage: lima togmx [OPTIONS]

Description:
    Build CHARMM27 coordinates, topology, and position restraints from a
    protein PDB or mmCIF file. Files are written beside the input structure.

Options:
    -f, -structure [path]
        Input .pdb or .cif file. <Required>

    -name, -n [string]
        Optional output basename. Without it, writes conf.gro, topol.top,
        and posre.itp. With it, writes NAME.gro, NAME.top, and NAME_posre.itp.

    -water, -water-model [tip3p|tip4p|tips3p|tip5p|spc|spce]
        Water topology to include for later solvation. Default: tip3p.

Example:
    lima togmx -f protein.cif -water tip3p
)";

    ArgParser parser(helpText);
    fs::path structurePath;
    std::string name;
    std::string water = "tip3p";
    parser.AddOption({ "-f", "-structure" }, true, structurePath);
    parser.AddOption({ "-name", "-n" }, false, name);
    parser.AddOption({ "-water", "-water-model" }, false, water);
    parser.Parse(argc, argv);

    if (!name.empty() && fs::path{ name }.has_parent_path()) {
        throw std::runtime_error("togmx output name must be a basename");
    }

    auto conversion = Programs::ToGmx(structurePath, Programs::ParseWaterModel(water));
    const fs::path directory = structurePath.parent_path().empty() ? fs::current_path() : structurePath.parent_path();
    const fs::path groPath = directory / (name.empty() ? "conf.gro" : name + ".gro");
    const fs::path topPath = directory / (name.empty() ? "topol.top" : name + ".top");

    const auto& molecules = conversion.topology.GetSystem().molecules;
    if (molecules.size() != conversion.positionRestraints.size()) {
        throw std::runtime_error("Converted molecule and position-restraint counts differ");
    }
    for (std::size_t i = 0; i < molecules.size(); ++i) {
        const std::string suffix = i == 0 ? std::string{} : "_" + molecules[i].name;
        const fs::path filename = name.empty()
            ? fs::path{ "posre" + suffix + ".itp" }
            : fs::path{ name + "_posre" + suffix + ".itp" };
        conversion.topology.moleculetypes.at(molecules[i].name)->positionRestraintsInclude = filename;
        conversion.positionRestraints[i].printToFile(directory / filename);
    }

    conversion.grofile.printToFile(groPath);
    conversion.topology.printToFile(topPath);
    return 0;
}
