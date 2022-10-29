#pragma once




#define LIMA_DEBUGMODE

//#define LIMA_SAFERUN		// Use this for?
//#define LIMA_VERBOSE


// --- These can be toggled to test the capabilities of LIMA without the two improperly implemented features --- //
#define ENABLE_SOLVENTS				// Enables Explicit Solvents
#define ENABLE_BLJV				// Ensures no LJ between bonded particles
//#define IGNORE_NEIGHBOR_OVERFLOW	// LJ list doesnt include neighbors when limit is reached.
// -------------------------------------------------------------------------------------------------------------- //




// Disable this for faster simulations. 
#define ENABLE_DISPLAY

// DIHEDRALS can no longer be disabled, as they are critical for BLJI

// Implicit Non-overlapping Work Distribution
#define INWD

// Debugging
const bool print_compound_positions = false;



// SImulation specifics

constexpr float BOX_LEN = 7.f;		// Must be > twice the len of largest compound
constexpr float BOX_LEN_HALF = BOX_LEN / 2.f;

const int N_SOLVATE_MOLECULES = 2000;	


const int NEIGHBORLIST_MAX_COMPOUNDS = 64;
//const int NEIGHBORLIST_MAX_SOLVENTS = 5120;
const int NEIGHBORLIST_MAX_SOLVENTS = 6144;
//const int SOLVENTS_BATCHSIZE = 64;







const int MAX_SOLVENTS = 0xFFFF;
const int SOLVENT_TESTLIMIT= MAX_SOLVENTS;
constexpr float CUTOFF = 0.8f;	//nm/
//const int MAX_ATOM_TYPES = 16;



// OPTIMIZATION PARAMETERS
const int MAX_COMPOUND_PARTICLES = 64;



const int THREADS_PER_SOLVENTBLOCK = 128;	
const int THREADS_PER_COMPOUNDBLOCK = MAX_COMPOUND_PARTICLES; // Must be >= max comp particles


const int SIMULATION_STEPS = 20000;

const int RAMPUP_STEPS = 1;


const int STEPS_PER_NLIST_UPDATE = 4;



const int STEPS_PER_RENDER = 1;



// THERMOSTAT PARAMETERS
const int STEPS_PER_LOGTRANSFER = 2;		// Must be >= 3
//const int STEPS_PER_TRAJTRANSFER = 100;
const int STEPS_PER_THERMOSTAT = 20;			// Must be >= 3
const bool ENABLE_BOXTEMP = true;
const bool APPLY_THERMOSTAT = true;										// Switch to using forcefield_host first
const bool PRINT_TEMP = false;


const int MAX_COMPOUNDS = 0xFF;
const int MAX_ATOMS = 1'000'000;
const int MAX_ATOMS_IN_RESIDUE = 32;		// TODO SET YP AGAIN






// Logging constants
const int N_DATAGAN_VALUES = 3;
const int STEPS_PER_TRAINDATATRANSFER = 100;








// Related to compound bridges
const int COMPOUNDBRIDGES_IN_BUNDLE = 96;
const int MAX_PARTICLES_IN_BRIDGE = 64;
const int MAX_SINGLEBONDS_IN_BRIDGE = 32;
const int MAX_ANGLEBONDS_IN_BRIDGE = 32;
const int MAX_DIHEDRALBONDS_IN_BRIDGE = 32;




