#pragma once



#define LIMA_DEBUGMODE
//#define LIMA_SAFERUN		// Use this for?
//#define LIMA_VERBOSE



// -------------------------------------------- Physics Parameters ---------------------------------------------- //
const int RAMPUP_STEPS = 0;				// Set to 0 to disable
constexpr float RAMPUP_MOMENTUM_SCALAR = 0.2f;
constexpr float MAX_RAMPUP_DIST = 0.0001f;	// [nm] how far any particle is max allowed to move during ramp-up

constexpr float VEL_RMS_SCALAR = 0.f;		// Set to 0 to freeze solvents
constexpr float CUTOFF = 0.9f;	//nm/
// -------------------------------------------------------------------------------------------------------------- //






// ------------------------------------------------ Box Parameters ---------------------------------------------- //
constexpr float BOX_LEN = 7.f;		// Must be > twice the len of largest compound
constexpr float BOX_LEN_HALF = BOX_LEN / 2.f;
// -------------------------------------------------------------------------------------------------------------- //






// -------------------------------------------- Simulation Parameters ------------------------------------------- //
//const int SIMULATION_STEPS = 180;
const bool print_compound_positions = false;		// what is tihs?
const bool DUMP_TRAJ = false;
const bool DUMP_POTE = false;
const bool POSTSIM_ANAL = true;
// -------------------------------------------------------------------------------------------------------------- //






// -------------------------------------------- Solvation Parameters -------------------------------------------- //
#define ENABLE_SOLVENTS				// Enables Explicit Solvents
const int MAX_SOLVENTS = 0xFFFF;
const int SOLVENT_TESTLIMIT = MAX_SOLVENTS;
const int N_SOLVATE_MOLECULES = 12000;			// Used when not loading from .conf file
// -------------------------------------------------------------------------------------------------------------- //






// ------------------------------------------ Optimization Parameters ------------------------------------------- //
const int MAX_COMPOUND_PARTICLES = 64;
const int MAX_COMPOUNDS = 0xFF;
const int MAX_ATOMS = 1'000'000;
const int MAX_ATOMS_IN_RESIDUE = 32;		// TODO SET YP AGAIN

const int NEIGHBORLIST_MAX_COMPOUNDS = 64;
const int NEIGHBORLIST_MAX_SOLVENTS = 6144;


// Related to compound bridges
const int COMPOUNDBRIDGES_IN_BUNDLE = 96;
const int MAX_PARTICLES_IN_BRIDGE = 64;
const int MAX_SINGLEBONDS_IN_BRIDGE = 32;
const int MAX_ANGLEBONDS_IN_BRIDGE = 32;
const int MAX_DIHEDRALBONDS_IN_BRIDGE = 32;
// -------------------------------------------------------------------------------------------------------------- //






// -------------------------------------------- Kernel Parameters ----------------------------------------------- //
const int THREADS_PER_SOLVENTBLOCK = 128;
const int THREADS_PER_COMPOUNDBLOCK = MAX_COMPOUND_PARTICLES;
// -------------------------------------------------------------------------------------------------------------- //





// ------------------------------------------- Temperature Parameters ------------------------------------------- //
const bool ENABLE_BOXTEMP = false;		// Calc box-temp
const bool APPLY_THERMOSTAT = false;		// Apply scalar based on temp	TODO: Switch to using forcefield_host first
const bool PRINT_TEMP = false;			// Force always print temp
const int STEPS_PER_THERMOSTAT = 10;			// Must be >= 3 why?
const int FIRST_TEMPERATURE_PRINT_STEP = RAMPUP_STEPS;
const int FIRST_THERMOSTAT_APPLICATION_STEP = RAMPUP_STEPS + 200;
constexpr float MAX_THERMOSTAT_SCALER = 0.01f / static_cast<float>(STEPS_PER_THERMOSTAT);
// -------------------------------------------------------------------------------------------------------------- //



// ------------------------------------------------ Display Parameters ---------------------------------------------- //
#define ENABLE_DISPLAY		// Disable this for faster simulations. 
const int STEPS_PER_RENDER = 5;
constexpr float FORCED_INTERRENDER_TIME = 0.f;		// [ms] Set to 0 for full speed sim
const int FIRST_INTERRENDER_WAIT = RAMPUP_STEPS;
// -------------------------------------------------------------------------------------------------------------- //


const int STEPS_PER_NLIST_UPDATE = 10;







// THERMOSTAT PARAMETERS
const int STEPS_PER_LOGTRANSFER = 5;		// Must be >= 3	why?
//const int STEPS_PER_TRAJTRANSFER = 100;

// Logging constants
const int N_DATAGAN_VALUES = 3;
const int STEPS_PER_TRAINDATATRANSFER = 100;













