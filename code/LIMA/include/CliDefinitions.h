#pragma once

#include <array>
#include <string_view>

#include "CliContext.h"

namespace Cli {

using CommandHandler = int(*)(const Context& ctx, int argc, char** argv);

// Handler-free so tools other than lima (eg. limaclitest) can enumerate the commands.
// lima.cpp maps each name to its handler.
struct CommandDefinition {
    std::string_view name;
    std::string_view summary;
    std::string_view helpText;
};

int RunMdrun(const Context& ctx, int argc, char** argv);
int RunBuildMembrane(const Context& ctx, int argc, char** argv);
int RunMakeSimParams(const Context& ctx, int argc, char** argv);
int RunSelfTest(const Context& ctx, int argc, char** argv);
int RunRender(const Context& ctx, int argc, char** argv);
int RunMakeBox(const Context& ctx, int argc, char** argv);
int RunSolvate(const Context& ctx, int argc, char** argv);
int RunInsertMolecule(const Context& ctx, int argc, char** argv);
int RunInsertMolecules(const Context& ctx, int argc, char** argv);
int RunEditConf(const Context& ctx, int argc, char** argv);
int RunEnergyMinimization(const Context& ctx, int argc, char** argv);
int RunToGmx(const Context& ctx, int argc, char** argv);
int RunDaemon(const Context& ctx, int argc, char** argv);
int RunWait(const Context& ctx, int argc, char** argv);

inline constexpr std::string_view MdrunHelp = R"(Usage: lima mdrun [OPTION]...

Run a molecular dynamics simulation.

Options:
  -c, --conf PATH              Input coordinates (default: ./conf.gro)
  -t, --topology PATH          Input topology (default: ./topol.top)
  -s, --simparams PATH         Simulation parameters (default: ./sim_params.txt)
      --set KEY=VALUE...       Override parameters of the simparams file, eg. --set ref_t=310
      --conf-out PATH          Output coordinates (default: ./out.gro)
      --trajectory PATH        Output trajectory, .trr or .xtc (compressed), written
                               during the run: a frame every trajectory_interval
                               steps (default 1000) and the final one. Disabled by
                               default
  -d, --display                Display the running simulation
      --async                  Run on the lima daemon and return at once; see 'lima wait'
      --uff                    Export trajectory.uff for ML training, every
                               data_logging_interval steps, kept in RAM for the run
  -h, --help                   Display this help and exit
)";

inline constexpr std::string_view BuildMembraneHelp = R"(Usage: lima buildmembrane [OPTION]...

Build a membrane with a specified lipid composition. Stockholm lipids 2020 are
supported by default; custom lipid files may be placed in the working directory.
Choose exactly one geometry option: --plane, --sphere, or --ellipsoid. If none
is given, LIMA creates a planar membrane at half the box height.

Options:
      --lipids NAME PERCENT... Lipid names and percentages (required)
  -b, --box-size NM [NM NM]    Cubic size or x/y/z dimensions (required)
      --plane Z                Planar bilayer mid-plane (default: half box height)
      --sphere X Y Z R         Vesicle center and bilayer mid-surface radius
      --ellipsoid X Y Z RX RY RZ
                                Ellipsoid center and three mid-surface radii
      --em-tolerance VALUE     EM force tolerance (default: 100 kJ/mol/nm)
      --seed INTEGER           Reproducible membrane construction seed (default: 0)
      --working-dir PATH       Output directory (default: current directory)
  -d, --display                Display the simulation
      --async                  Run on the lima daemon and return at once; see 'lima wait'
  -h, --help                   Display this help and exit

Example:
  lima buildmembrane --lipids DPPC 60 DOPC 40 --box-size 10
  lima buildmembrane --lipids DMPC 100 --box-size 16 --sphere 8 8 8 5
  lima buildmembrane --lipids POPC 70 cholesterol 30 --box-size 24 \
    --ellipsoid 12 12 12 7 8 9 --seed 42
  lima buildmembrane --lipids DOPC 100 --box-size 12 --plane 5.5
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
  -x, --trajectory PATH        Play back an .xtc or .trr trajectory of the structure
                               (space plays/pauses, comma and period step a frame)
  -r, --representation NAME    Initial representation: atomname, charge, atomid,
                               compound, backbone, lipids, licorice or surface
                               (default: atomname)
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

inline constexpr std::string_view SolvateHelp = R"(Usage: lima solvate [OPTION]...

Fill the unoccupied volume of a simulation box with water and update its
topology.

Options:
  -c, --conf PATH              Input coordinates (default: ./conf.gro)
  -t, --topology PATH          Input topology (default: ./topol.top)
      --conf-out PATH          Output coordinates (default: INPUT_solvated.gro)
      --topology-out PATH      Output topology (default: INPUT_solvated.top)
  -p, --density COUNT          Target waters per nm^3 (default: 34)
      --water-model MODEL      Water model name or .itp file (default: spce)
  -h, --help                   Display this help and exit

A model name uses that model's .itp from the system's forcefield, or from
CHARMM27. An .itp file must define one 3-atom water moleculetype, with
[ settles ] or with bonds and an angle that have explicit parameters.
Water models: tip3p, tips3p, spc, spce (tip4p and tip5p need virtual sites,
which LIMA does not support)

Example:
  lima solvate --conf membrane.gro --topology membrane.top
  lima solvate --conf protein.gro --topology protein.top --water-model tip3p
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

Transform a single molecule, or set the box of a whole system. Operations are
applied in the order shown below.

Options:
  -c, --conf PATH              Input coordinates (required)
  -t, --topology PATH          Input topology (required; not modified)
      --conf-out PATH          Output coordinates (default: overwrite input)
  -w, --whole                  Make a PBC-fragmented molecule whole
      --set-center X Y Z       Set the geometric center
  -r, --rotate X Y Z           Rotate around z, y, then x, in radians
  -b, --box-size NM [NM NM]    Set a cubic or x/y/z box and center the system in it
  -d, --distance NM            Set a cubic box of the system's diameter (largest
                               atom-atom distance) plus twice NM, and center the
                               system in it, as GROMACS' editconf -d -bt cubic
  -h, --help                   Display this help and exit

--whole, --set-center and --rotate need a single molecule; --box-size and
--distance apply to every atom, eg. a multi-chain protein from togmx.

Example:
  lima togmx -f 1ubq.pdb
  lima editconf -c conf.gro -t topol.top -d 1.0
  lima solvate -c conf.gro -t topol.top
)";

inline constexpr std::string_view EnergyMinimizationHelp = R"(Usage: lima em [OPTION]...

Energy-minimize a simulation using standard parameters.

Options:
  -c, --conf PATH              Input coordinates (default: ./conf.gro)
  -t, --topology PATH          Input topology (default: ./topol.top)
      --conf-out PATH          Output coordinates (default: ./conf.gro)
      --em-tolerance VALUE     Force tolerance (default: 100 kJ/mol/nm)
  -d, --display                Display the simulation
      --async                  Run on the lima daemon and return at once; see 'lima wait'
  -h, --help                   Display this help and exit
)";

inline constexpr std::string_view DaemonHelp = R"(Usage: lima daemon status|stop

GPU commands (mdrun, em and buildmembrane) started with --async run in a
background lima process, the daemon, so commands submitted together share the
GPU efficiently. Without --async a command runs in its own process:

  for d in run1 run2 run3; do lima mdrun --async -c $d/conf.gro -t $d/topol.top \
    --conf-out $d/out.gro; done
  lima wait

The daemon starts by itself and exits after 2 minutes without work. With -d,
the simulations it runs together are shown side by side in one window.

Commands:
  status                       List the commands the daemon is running
  stop                         Let running commands finish, then exit

Environment:
  LIMA_NO_DAEMON=1             Disable the daemon (--async then fails)
  LIMA_DAEMON_IDLE=SECONDS     Idle time before the daemon exits (default: 120)
  LIMA_DAEMON_DIR=PATH         Where the daemon keeps its socket and log
)";

inline constexpr std::string_view WaitHelp = R"(Usage: lima wait

Wait until every command started with --async has finished, then list each with
its result. Exits 1 if any of them failed. The output of a failed command is kept
in lima-JOB.log in the directory it was started from; successful ones leave none.

Example:
  for d in run1 run2 run3; do lima mdrun --async -c $d/conf.gro -t $d/topol.top \
    --conf-out $d/out.gro; done
  lima wait

Options:
  -h, --help                   Display this help and exit
)";

inline constexpr std::string_view ToGmxHelp = R"(Usage: lima togmx [OPTION]...

Build CHARMM27 or AMBER99SB coordinates, topology, and position restraints from a
protein PDB or mmCIF file. Prints what it decided (disulfides, histidine states) and
what it skipped (HETATM records such as waters, ions and ligands).

Options:
  -f, --structure PATH         Input .pdb or .cif file (required)
  -o, --output-dir DIR         Where to write the files (default: beside the input)
  -n, --name NAME              Optional output basename
      --ff FORCEFIELD          charmm27 or amber99sb (default: charmm27). AMBER99SB only
                               has charged termini
      --water-model MODEL      Water topology (default: tip3p)
      --chains A,B             Only convert these chains (default: all)
      --altloc ID              Alternate location to keep (default: A)
      --his MODE               hsd, hse, hsp or auto (default: hsd). Input names
                               HISD, HISE, HISH always win
      --ter-n MODE             charged (NH3+) or neutral (NH2) (default: charged)
      --ter-c MODE             charged (COO-), neutral (COOH), amide (CONH2),
                               or methylamide (CONHCH3) (default: charged)
      --no-disulfides          Don't bond cysteines whose SG atoms are 0.17-0.24 nm apart
      --ignore-h               Ignore the input's hydrogens and build them all
  -h, --help                   Display this help and exit

Water models: tip3p, tip4p, tips3p (charmm27 only), tip5p, spc, spce
Disulfides between chains are not bonded (each chain is its own molecule); they
are reported as a warning.
)";

inline constexpr std::array Commands{
    CommandDefinition{ "mdrun", "Run a molecular dynamics simulation.", MdrunHelp },
    CommandDefinition{ "buildmembrane", "Build a membrane structure.", BuildMembraneHelp },
    CommandDefinition{ "makesimparams", "Write default simulation parameters.", MakeSimParamsHelp },
    CommandDefinition{ "selftest", "Run LIMA's internal self-test.", SelfTestHelp },
    CommandDefinition{ "render", "Render a molecular structure.", RenderHelp },
    CommandDefinition{ "makebox", "Create an empty simulation box.", MakeBoxHelp },
    CommandDefinition{ "solvate", "Fill a simulation box with water.", SolvateHelp },
    CommandDefinition{ "insertmolecule", "Insert one molecule into a box.", InsertMoleculeHelp },
    CommandDefinition{ "insertmolecules", "Insert multiple molecules into a box.", InsertMoleculesHelp },
    CommandDefinition{ "editconf", "Transform molecular coordinates.", EditConfHelp },
    CommandDefinition{ "em", "Energy-minimize a simulation.", EnergyMinimizationHelp },
    CommandDefinition{ "togmx", "Convert a protein structure to CHARMM27 or AMBER99SB files.", ToGmxHelp },
    CommandDefinition{ "daemon", "Show or stop the background GPU process.", DaemonHelp },
    CommandDefinition{ "wait", "Wait for commands started with --async.", WaitHelp },
};

} // namespace Cli
