#pragma once

#include <math.h>
#include <cstdint>
#include <limits.h>


// LIMASAFEMODE slightly alters the outcome of sims. Even overwrite enabling it in impropers, for a
// sim with no impropers has this effect. It is very wierd, and i fear i have some undefined behavior
// somewhere in the code
//#define LIMASAFEMODE
//#define LIMAPUSH
#if defined LIMAPUSH && defined LIMASAFEMODE
#error These are mutually exclusive
#endif

//#define FORCE_NAN_CHECK
const bool LIMA_PUSH = true;

#define ENABLE_LJ
#define ENABLE_INTEGRATEPOSITION

const bool ENABLE_ES_SR = true;
const bool ENABLE_ES_LR = false; // This is not deterministic, due to the atomicAdd in DistributeChargesToChargegrid

//#define GENERATETRAINDATA

//#define LIMAKERNELDEBUGMODE
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
constexpr int64_t NANO_TO_LIMA_i = static_cast<int64_t>(NANO_TO_LIMA);
constexpr float LIMA_TO_NANO = 1.f / NANO_TO_LIMA;
const int PICO_TO_LIMA = static_cast<int>(FEMTO_TO_LIMA) * 1000;

static_assert(NANO_TO_LIMA < INT_MAX/4, "LIMA Scale is so small it can create dangerous bugs");

constexpr float KILO = 1000.f;
constexpr double GIGA = 1e9;
constexpr double NANO = 1e-9;
constexpr double FEMTO = 1e-15;
constexpr double LIMA = NANO / NANO_TO_LIMA;

constexpr float kcalToJoule = 4184.f;
constexpr float degreeToRad = 2.f * PI / 360.f;
constexpr float AngToNm = 0.1f;
const float rminToSigma = 1.f / powf(2.f, (1.f / 6.f));


const float DEG_TO_RAD = PI / 180.f;

const int MAX_REPRESENTABLE_DIFF_NM = 16;	// I should probably do this some other way..

constexpr double BOLTZMANNCONSTANT = 1.38066e-23f;	// [J/K]
constexpr double AVOGADROSNUMBER = 6.02214076e23;	
constexpr double COULOMBCONSTANT = 8.9875517873681764e9;	// [N m^2 / C^2]
constexpr double ELEMENTARYCHARGE = 1.602176634e-19;	// [C]

constexpr float elementaryChargeToKiloCoulombPerMole = ELEMENTARYCHARGE * AVOGADROSNUMBER / KILO;
// -------------------------------------------------------------------------------------------------------------- //





// -------------------------------------------- Solvation Parameters -------------------------------------------- //
#define ENABLE_SOLVENTS				// Enables Explicit Solvents
const size_t MAX_SOLVENTS = INT32_MAX-1;	// limited by boxparams
constexpr float DEFAULT_TINYMOL_START_TEMPERATURE = 310.f;	// [K]
// -------------------------------------------------------------------------------------------------------------- //






// ------------------------------------------ Optimization Parameters ------------------------------------------- //
const bool HARD_CUTOFF = true;
const bool ENABLE_POTE = true;
const bool IGNORE_HYDROGEN = false;
const int GRIDNODE_QUERY_RANGE = 2;

// If we go larger, a single compound can stretch over 2 nm!
//constexpr int MAX_COMPOUND_PARTICLES = IGNORE_HYDROGEN ? 48 : 64;
constexpr int MAX_COMPOUND_PARTICLES = 64;
const int MAX_COMPOUNDS = UINT16_MAX-1;			// Arbitrary i think. true max int16_t max - 1. Can also cause trouble when the bondedparticlesLUT static array becomes very large bytewise..

const int NEIGHBORLIST_MAX_COMPOUNDS = 512;	// TODO: We need to work on getting this number down!
//const int NEIGHBORLIST_MAX_SOLVENTS = 6144;

const bool USE_ATOMICS_FOR_BONDS_RESULTS = false;


// Related to compound bridges
const int MAX_COMPOUNDBRIDGES = MAX_COMPOUNDS;	// Wtf is this param?
const int MAX_PARTICLES_IN_BRIDGE = 32+16;	// Limited to 255 by getBridgelocalIdOfParticle, since id 255 is invalid
const int MAX_SINGLEBONDS_IN_BRIDGE = 8;
const int MAX_ANGLEBONDS_IN_BRIDGE = 16;
const int MAX_DIHEDRALBONDS_IN_BRIDGE = 64 + 16 + 16;
const int MAX_IMPROPERDIHEDRALBONDS_IN_BRIDGE = 4;
const int MAX_COMPOUNDS_IN_BRIDGE = 4;	// Some bridges span more than 2 compounds, for example the loop between beta plates

const int MAX_SAFE_SHIFT = 6;	// Maxmimum manhattan dist that it is safe to shift








// -------------------------------------------- Kernel Parameters ----------------------------------------------- //
const int THREADS_PER_COMPOUNDBLOCK = MAX_COMPOUND_PARTICLES;
// -------------------------------------------------------------------------------------------------------------- //


// -------------------------------------------- Neighborlist Parameters ----------------------------------------- //
const int STEPS_PER_NLIST_UPDATE = 5;
const bool ALLOW_ASYNC_NLISTUPDATE = true;
// -------------------------------------------------------------------------------------------------------------- //
