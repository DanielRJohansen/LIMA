#pragma once

#include "Programs.h"
#include "argparser.h"

#include <stdexcept>

inline int cif2gmx(int argc, char** argv) {
    const std::string helpText = R"(
Usage: lima cif2gmx [OPTIONS]

Description:
    Build CHARMM27 coordinates, topology, and position restraints from a
    single-chain protein mmCIF file. Files are written beside the input CIF.

Options:
    -f, -cif [path]
        Input mmCIF file. <Required>

    -name, -n [string]
        Optional output basename. Without it, writes conf.gro, topol.top,
        and posre.itp. With it, writes NAME.gro, NAME.top, and NAME_posre.itp.

    -water, -water-model [tip3p|tip4p|tips3p|tip5p|spc|spce]
        Water topology to include for later solvation. Default: tip3p.

    -help, -h
        Display this help text and exit.

Example:
    lima cif2gmx -f protein.cif -water tip3p
)";

    ArgParser parser(helpText);
    fs::path cifPath;
    std::string name;
    std::string water = "tip3p";
    parser.AddOption({ "-f", "-cif" }, true, cifPath);
    parser.AddOption({ "-name", "-n" }, false, name);
    parser.AddOption({ "-water", "-water-model" }, false, water);
    parser.Parse(argc, argv);

    Programs::WaterModel waterModel;
    if (water == "tip3p") waterModel = Programs::WaterModel::Tip3p;
    else if (water == "tip4p") waterModel = Programs::WaterModel::Tip4p;
    else if (water == "tips3p") waterModel = Programs::WaterModel::Tips3p;
    else if (water == "tip5p") waterModel = Programs::WaterModel::Tip5p;
    else if (water == "spc") waterModel = Programs::WaterModel::Spc;
    else if (water == "spce") waterModel = Programs::WaterModel::Spce;
    else throw std::runtime_error("Unsupported CHARMM27 water model: " + water);

    Programs::cif2gmx(cifPath, name.empty() ? std::nullopt : std::optional<std::string>{ name }, waterModel);
    return 0;
}
