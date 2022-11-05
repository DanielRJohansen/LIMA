#pragma once



#define LIMA_DEBUGMODE
//#define LIMA_SAFERUN		// Use this for?
//#define LIMA_VERBOSE



// -------------------------------------------- Physics Parameters ---------------------------------------------- //
const bool INTEGRATION_RAMPUP_ENABLED = 1;
const int RAMPUP_STEPS = 500;
constexpr float RAMPUP_MOMENTUM_SCALAR = 0.2f;

constexpr float VEL_RMS_SCALAR = 1.f;		// Set to 0 to freeze solvents
constexpr float CUTOFF = 0.9f;	//nm/
// -------------------------------------------------------------------------------------------------------------- //






// ------------------------------------------------ Box Parameters ---------------------------------------------- //
constexpr float BOX_LEN = 7.f;		// Must be > twice the len of largest compound
constexpr float BOX_LEN_HALF = BOX_LEN / 2.f;
constexpr float FORCED_INTERRENDER_TIME = 0.f;		// [ms] Set to 0 for full speed sim
// -------------------------------------------------------------------------------------------------------------- //






// -------------------------------------------- Simulation Parameters ------------------------------------------- //
const int SIMULATION_STEPS = 20000;
const bool print_compound_positions = false;		// what is tihs?
// -------------------------------------------------------------------------------------------------------------- //






// -------------------------------------------- Solvation Parameters -------------------------------------------- //
#define ENABLE_SOLVENTS				// Enables Explicit Solvents
const int MAX_SOLVENTS = 0xFFFF;
const int SOLVENT_TESTLIMIT = MAX_SOLVENTS;
const int N_SOLVATE_MOLECULES = 8000;			// Used when not loading from .conf file
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
const bool ENABLE_BOXTEMP = true;		// Calc box-temp
const bool APPLY_THERMOSTAT = true;		// Apply scalar based on temp	TODO: Switch to using forcefield_host first
const bool PRINT_TEMP = false;			// Force always print temp
const int STEPS_PER_THERMOSTAT = 5;			// Must be >= 3 why?
const int FIRST_TEMPERATURE_PRINT_STEP = 0;// RAMPUP_STEPS* INTEGRATION_RAMPUP_ENABLED;
const int FIRST_THERMOSTAT_APPLICATION_STEP = 200;
constexpr float MAX_THERMOSTAT_SCALER = 0.001f;
// -------------------------------------------------------------------------------------------------------------- //



// ------------------------------------------------ Display Parameters ---------------------------------------------- //
#define ENABLE_DISPLAY		// Disable this for faster simulations. 
const int STEPS_PER_RENDER = 80;
// -------------------------------------------------------------------------------------------------------------- //


const int STEPS_PER_NLIST_UPDATE = 10;







// THERMOSTAT PARAMETERS
const int STEPS_PER_LOGTRANSFER = 5;		// Must be >= 3	why?
//const int STEPS_PER_TRAJTRANSFER = 100;

// Logging constants
const int N_DATAGAN_VALUES = 3;
const int STEPS_PER_TRAINDATATRANSFER = 100;













