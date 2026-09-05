#pragma once

#include <array>
#include <string_view>

namespace Cli {

using CommandHandler = int(*)(int argc, char** argv);

struct CommandDefinition {
    std::string_view name;
    std::string_view summary;
    std::string_view helpText;
    CommandHandler handler;
};

int RunMdrun(int argc, char** argv);
int RunBuildMembrane(int argc, char** argv);
int RunMakeSimParams(int argc, char** argv);
int RunSelfTest(int argc, char** argv);
int RunRender(int argc, char** argv);
int RunMakeBox(int argc, char** argv);
int RunInsertMolecule(int argc, char** argv);
int RunInsertMolecules(int argc, char** argv);
int RunEditConf(int argc, char** argv);
int RunEnergyMinimization(int argc, char** argv);
int RunToGmx(int argc, char** argv);

inline constexpr std::string_view MdrunHelp = R"(Usage: lima mdrun [OPTION]...

Run a molecular dynamics simulation.

Options:
  -c, --conf PATH              Input coordinates (default: ./conf.gro)
  -t, --topology PATH          Input topology (default: ./topol.top)
  -s, --simparams PATH         Simulation parameters (default: ./sim_params.txt)
      --conf-out PATH          Output coordinates (default: ./out.gro)
      --trajectory PATH        Output trajectory; disabled by default
  -d, --display                Display the running simulation
      --uff                    Export trajectory.uff for ML training
  -h, --help                   Display this help and exit
)";

inline constexpr std::string_view BuildMembraneHelp = R"(Usage: lima buildmembrane [OPTION]...

Build a membrane with a specified lipid composition. Stockholm lipids 2020 are
supported by default; custom lipid files may be placed in the working directory.

Options:
      --lipids NAME PERCENT... Lipid names and percentages (required)
  -b, --box-size NM [NM NM]    Cubic size or x/y/z dimensions (required)
  -c, --center-z NM            Membrane center (default: half the box height)
      --em-tolerance VALUE     EM force tolerance (default: 100 kJ/mol/nm)
      --working-dir PATH       Output directory (default: current directory)
  -d, --display                Display the simulation
  -h, --help                   Display this help and exit

Example:
  lima buildmembrane --lipids DPPC 60 DOPC 40 --box-size 10
)";

inline constexpr std::string_view MakeSimParamsHelp = R"(Usage: lima makesimparams [OPTION]...

Write default simulation parameters to sim_params.txt in the current directory.

Options:
  -h, --help                   Display this help and exit
)";

inline constexpr std::string_view SelfTestHelp = R"(Usage: lima selftest [OPTION]...

Run LIMA's internal self-test.

Options:
  -h, --help                   Display this help and exit
)";

inline constexpr std::string_view RenderHelp = R"(Usage: lima render [OPTION]...

Render a .gro file, or convert and render a protein .pdb/.cif file.

Options:
  -f, --structure PATH         Input structure (default: ./conf.gro)
  -t, --topology PATH          Topology for .gro input (default: ./topol.top)
      --water-model MODEL      Conversion water model (default: tip3p)
  -w, --whole                  Unwrap and center the molecule
      --hide-water             Hide water molecules
      --highlight INDEX...     Highlight zero-based atom indices
  -h, --help                   Display this help and exit

Water models: tip3p, tip4p, tips3p, tip5p, spc, spce
)";

inline constexpr std::string_view MakeBoxHelp = R"(Usage: lima makebox [OPTION]...

Create an empty cubic coordinate box and topology.

Options:
  -b, --box-size NM            Cubic box size (required)
  -n, --name NAME              Output basename (default: conf/topol)
  -h, --help                   Display this help and exit

Example:
  lima makebox --name mybox --box-size 14
)";

inline constexpr std::string_view InsertMoleculeHelp = R"(Usage: lima insertmolecule [OPTION]...

Insert one molecule into a target box.

Options:
      --conf-source PATH       Source coordinates (required)
      --top-source PATH        Source topology (required)
      --conf-target PATH       Target coordinates (default: ./conf.gro)
      --top-target PATH        Target topology (default: ./topol.top)
  -p, --position X Y Z         Insertion position (default: box center)
  -h, --help                   Display this help and exit
)";

inline constexpr std::string_view InsertMoleculesHelp = R"(Usage: lima insertmolecules [OPTION]...

Insert multiple uniformly distributed molecules into a target box.

Options:
      --conf-source PATH       Source coordinates (required)
      --top-source PATH        Source topology (required)
      --conf-target PATH       Target coordinates (default: ./conf.gro)
      --top-target PATH        Target topology (default: ./topol.top)
  -n, --num-insertions COUNT   Number of insertions (required)
      --rotate-randomly        Randomly rotate each molecule
  -d, --display                Display energy minimization
  -h, --help                   Display this help and exit
)";

inline constexpr std::string_view EditConfHelp = R"(Usage: lima editconf [OPTION]...

Transform a single molecule. Operations are applied in the order shown below.

Options:
  -c, --conf PATH              Input coordinates (required)
  -t, --topology PATH          Input topology (required; not modified)
      --conf-out PATH          Output coordinates (default: overwrite input)
  -w, --whole                  Make a PBC-fragmented molecule whole
      --set-center X Y Z       Set the geometric center
  -r, --rotate X Y Z           Rotate around z, y, then x, in radians
  -h, --help                   Display this help and exit
)";

inline constexpr std::string_view EnergyMinimizationHelp = R"(Usage: lima em [OPTION]...

Energy-minimize a simulation using standard parameters.

Options:
  -c, --conf PATH              Input coordinates (default: ./conf.gro)
  -t, --topology PATH          Input topology (default: ./topol.top)
      --conf-out PATH          Output coordinates (default: ./conf.gro)
      --em-tolerance VALUE     Force tolerance (default: 100 kJ/mol/nm)
  -d, --display                Display the simulation
  -h, --help                   Display this help and exit
)";

inline constexpr std::string_view ToGmxHelp = R"(Usage: lima togmx [OPTION]...

Build CHARMM27 coordinates, topology, and position restraints from a protein
PDB or mmCIF file. Files are written beside the input structure.

Options:
  -f, --structure PATH         Input .pdb or .cif file (required)
  -n, --name NAME              Optional output basename
      --water-model MODEL      Water topology (default: tip3p)
  -h, --help                   Display this help and exit

Water models: tip3p, tip4p, tips3p, tip5p, spc, spce
)";

inline constexpr std::array Commands{
    CommandDefinition{ "mdrun", "Run a molecular dynamics simulation.", MdrunHelp, RunMdrun },
    CommandDefinition{ "buildmembrane", "Build a membrane structure.", BuildMembraneHelp, RunBuildMembrane },
    CommandDefinition{ "makesimparams", "Write default simulation parameters.", MakeSimParamsHelp, RunMakeSimParams },
    CommandDefinition{ "selftest", "Run LIMA's internal self-test.", SelfTestHelp, RunSelfTest },
    CommandDefinition{ "render", "Render a molecular structure.", RenderHelp, RunRender },
    CommandDefinition{ "makebox", "Create an empty simulation box.", MakeBoxHelp, RunMakeBox },
    CommandDefinition{ "insertmolecule", "Insert one molecule into a box.", InsertMoleculeHelp, RunInsertMolecule },
    CommandDefinition{ "insertmolecules", "Insert multiple molecules into a box.", InsertMoleculesHelp, RunInsertMolecules },
    CommandDefinition{ "editconf", "Transform molecular coordinates.", EditConfHelp, RunEditConf },
    CommandDefinition{ "em", "Energy-minimize a simulation.", EnergyMinimizationHelp, RunEnergyMinimization },
    CommandDefinition{ "togmx", "Convert a protein structure to CHARMM27 files.", ToGmxHelp, RunToGmx },
};

} // namespace Cli
