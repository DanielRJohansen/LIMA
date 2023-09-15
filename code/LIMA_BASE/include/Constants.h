#pragma once

#include <math.h>
#include <cstdint>
#include <limits.h>

#include "UserConstants.h"

#define LIMASAFEMODE
//#define LIMAPUSH
#if defined LIMAPUSH && defined LIMASAFEMODE
#error These are mutually exclusive
#endif


#define LIMAKERNELDEBUGMODE
//#define DONTGENDATA

constexpr float PI = 3.14159f;


// -------------------------------------------- Physics Parameters ---------------------------------------------- //
const int RAMPUP_STEPS = 0;					// Set to 0 to disable
constexpr float RAMPUP_MOMENTUM_SCALAR = 0.2f;
constexpr float MAX_RAMPUP_DIST = 0.0001f;	// [nm] how far any particle is max allowed to move during ramp-up

constexpr float VEL_RMS_SCALAR = 0.f;		// Set to 0 to freeze solvents

//constexpr float LIMA_SCALE = 1.f;// 1e-6f;			// size of 1 lima unit in nm or ns or whatever
constexpr float NANO_TO_FEMTO = 1e+6f;				// Allow for quickly changing all units from femto to another
constexpr float NANO_TO_PICO = 1e+3f;
constexpr float FEMTO_TO_LIMA = 100.f;		// >>7 to get fm when uint
constexpr float LIMA_TO_FEMTO = 1.f / FEMTO_TO_LIMA;

constexpr float NANO_TO_LIMA = FEMTO_TO_LIMA * NANO_TO_FEMTO;
constexpr float LIMA_TO_NANO = 1.f / NANO_TO_LIMA;
const int PICO_TO_LIMA = static_cast<int>(FEMTO_TO_LIMA) * 1000;

static_assert(NANO_TO_LIMA < INT_MAX/4, "LIMA Scale is so small it can create dangerous bugs");

constexpr float kcalToJoule = 4184.f;
constexpr float degreeToRad = 2.f * PI / 360.f;
constexpr float AngToNm = 0.1f;
const float rminToSigma = powf(2.f, (1.f / 6.f));

const int MAX_REPRESENTABLE_DIFF_NM = 16;	// I should probably do this some other way..

constexpr float CUTOFF_LM = CUTOFF_NM * NANO_TO_LIMA;				// fm

constexpr double BOLTZMANNCONSTANT = 1.38066e-23f;	// [J/K]
constexpr double AVOGADROSNUMBER = 6.02214076e23;	
// -------------------------------------------------------------------------------------------------------------- //






// ------------------------------------------------ Box Parameters ---------------------------------------------- //
//constexpr int _BOX_LEN_PM = 18000;
constexpr float BOX_LEN_NM = static_cast<float>(_BOX_LEN_PM) / 1000.f;

const int64_t BOX_LEN_i = static_cast<std::int64_t>(_BOX_LEN_PM) * PICO_TO_LIMA;
constexpr float BOX_LEN = BOX_LEN_NM * NANO_TO_LIMA;		// Must be > twice the len of largest compound
constexpr float BOX_LEN_HALF = BOX_LEN / 2.f;
constexpr float BOX_LEN_HALF_NM = BOX_LEN_NM / 2.f;

constexpr float BOX_LEN_FM = BOX_LEN * LIMA_TO_FEMTO;
constexpr float BOX_LEN_HALF_FM = BOX_LEN_FM / 2.f;

// -------------------------------------------------------------------------------------------------------------- //






// -------------------------------------------- Simulation Parameters ------------------------------------------- //
const bool print_compound_positions = false;		// what is tihs?
const bool DUMP_TRAJ = false;		// Make sure this is correct with DONTGETDATA flag
const bool DUMP_POTE = true;
const bool POSTSIM_ANAL = true;
// -------------------------------------------------------------------------------------------------------------- //






// -------------------------------------------- Solvation Parameters -------------------------------------------- //
#define ENABLE_SOLVENTS				// Enables Explicit Solvents
const size_t MAX_SOLVENTS = INT32_MAX-1;	// limited by boxparams

const int MAX_SOLVENTS_IN_BLOCK = 256;
const int STEPS_PER_SOLVENTBLOCKTRANSFER = 5;	// If we go below 2, we might see issue in solventtransfers
const int SOLVENTBLOCK_TRANSFERSTEP = STEPS_PER_SOLVENTBLOCKTRANSFER - 1;
// -------------------------------------------------------------------------------------------------------------- //






// ------------------------------------------ Optimization Parameters ------------------------------------------- //
const bool HARD_CUTOFF = true;
const bool CALC_POTE = true;
const bool IGNORE_HYDROGEN = false;
const int GRIDNODE_QUERY_RANGE = 2;


const int MAX_COMPOUND_PARTICLES = 48;	// If we go larger, a single compound can stretch over 2 nm!
const int MAX_COMPOUNDS = 2048;			// Arbitrary i think. true max int16_t max - 1. Can also cause trouble when the bondedparticlesLUT static array becomes very large bytewise..

const int NEIGHBORLIST_MAX_COMPOUNDS = 256;
const int NEIGHBORLIST_MAX_SOLVENTS = 6144;


// Related to compound bridges
const int MAX_COMPOUNDBRIDGES = MAX_COMPOUNDS;	// Wtf is this param?
const int MAX_PARTICLES_IN_BRIDGE = 32;
const int MAX_SINGLEBONDS_IN_BRIDGE = 4;
const int MAX_ANGLEBONDS_IN_BRIDGE = 16;
const int MAX_DIHEDRALBONDS_IN_BRIDGE = 32;
const int MAX_IMPROPERDIHEDRALBONDS_IN_BRIDGE = 4;
const int MAX_COMPOUNDS_IN_BRIDGE = 4;	// Some bridges span more than 2 compounds, for example the loop between beta plates

const int MAX_SAFE_SHIFT = 6;	// Maxmimum manhattan dist that it is safe to shift

// Related to forcefield / constant memory
const int MAX_ATOM_TYPES = 48;

constexpr float MAX_COMPOUND_RADIUS = 1.9f;	// was 1.5
// -------------------------------------------------------------------------------------------------------------- //






// -------------------------------------------- Kernel Parameters ----------------------------------------------- //
const int THREADS_PER_COMPOUNDBLOCK = MAX_COMPOUND_PARTICLES;
// -------------------------------------------------------------------------------------------------------------- //





// ------------------------------------------- Temperature Parameters ------------------------------------------- //
const bool ENABLE_BOXTEMP	= true;		// Calc box-temp
const bool PRINT_TEMP = false;			// Force always print temp
const int STEPS_PER_THERMOSTAT = 200;			// Must be >= 3 why?
constexpr float MAX_THERMOSTAT_SCALER = 0.001f / static_cast<float>(STEPS_PER_THERMOSTAT);	// change vel by 0.1% over NSTEPS
// -------------------------------------------------------------------------------------------------------------- //



// ------------------------------------------------ Display Parameters ---------------------------------------------- //
#define ENABLE_DISPLAY		// Disable this for faster simulations. 
#if defined ENABLE_DISPLAY && defined __LINUX__
#error It is not allowed to use display on linux as of right now
#endif

const int STEPS_PER_RENDER = 50;
constexpr float FORCED_INTERRENDER_TIME = 0.f;		// [ms] Set to 0 for full speed sim
// -------------------------------------------------------------------------------------------------------------- //

// -------------------------------------------- Neighborlist Parameters ----------------------------------------- //
const int STEPS_PER_NLIST_UPDATE = 10;
const bool ALLOW_ASYNC_NLISTUPDATE = true;
// -------------------------------------------------------------------------------------------------------------- //






const int LOG_EVERY_N_STEPS = 5;
const int STEPS_PER_LOGTRANSFER =  10;		// Must be >= 3	why?
static_assert(STEPS_PER_LOGTRANSFER% LOG_EVERY_N_STEPS == 0, "Log intervals doesn't match");

// Logging constants
const int N_DATAGAN_VALUES = 3;
const int STEPS_PER_TRAINDATATRANSFER = 100;
