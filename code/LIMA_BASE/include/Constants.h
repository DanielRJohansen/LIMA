#pragma once

#include <cstdint>
#include <math.h>


// -------------------------------------------- Debug Parameters -------------------------------------------- //
constexpr bool INDEXING_CHECKS = false;
constexpr bool SYNC_ALL_KERNELS = false;	// Disallow async/concurrent kernels
constexpr bool FORCE_CHECKS = false;		// Check force is not NaN or Inf
constexpr bool POSITION_CHECKS = false;		// Check if near overflow when switching to int representation


//#define FORCE_NAN_CHECK

#define ENABLE_LJ
#define ENABLE_INTEGRATEPOSITION

const bool ENABLE_ES_SR = true;
const bool ENABLE_ES_LR = true; // This is not deterministic, due to the atomicAdd in DistributeChargesToChargegrid

const bool ENABLE_ERFC_FOR_EWALD = true;

const bool ENABLE_UREYBRADLEY = true;

const bool USE_PRECOMPUTED_BSPLINES = false;
const bool USE_PRECOMPUTED_ERFCSCALARS = false;
//#define GENERATETRAINDATA

//#define LIMAKERNELDEBUGMODE
//#define DONTGENDATA






// -------------------------------------------- Physics Parameters ---------------------------------------------- //
constexpr float KILO = 1000.f;
constexpr double GIGA = 1e9;
constexpr double NANO = 1e-9;
constexpr double FEMTO = 1e-15;

constexpr double NANO_TO_FEMTO = 1e6;
constexpr double FEMTO_TO_NANO = 1e-6;
constexpr double NANO_TO_PICO = 1e3;

constexpr float PI = 3.14159f;
constexpr float kcalToJoule = 4184.f;
constexpr float degreeToRad = 2.f * PI / 360.f;
constexpr float AngToNm = 0.1f;
const float rminToSigma = 1.f / powf(2.f, (1.f / 6.f));

const float DEG_TO_RAD = PI / 180.f;

const int MAX_REPRESENTABLE_DIFF_NM = 16;	// I should probably do this some other way..

constexpr double BOLTZMANNCONSTANT = 1.38066e-23f;	// [J/K]
constexpr double AVOGADROSNUMBER = 6.02214076e23;	
constexpr double COULOMBCONSTANT = 8.9875517873681764e9;	// [n*m^2/C^2] == [ J*m / C^2]
constexpr double ELEMENTARYCHARGE = 1.602176634e-19;	// [C]

constexpr float elementaryChargeToKiloCoulombPerMole = ELEMENTARYCHARGE * AVOGADROSNUMBER / KILO;
// -------------------------------------------------------------------------------------------------------------- //





// -------------------------------------------- Solvation Parameters -------------------------------------------- //
#define ENABLE_SOLVENTS				// Enables Explicit Solvents
const size_t MAX_SOLVENTS = INT32_MAX-1;	// limited by boxparams
constexpr float DEFAULT_TINYMOL_START_TEMPERATURE = 310.f;	// [K]
constexpr bool AllAtom = true;
// -------------------------------------------------------------------------------------------------------------- //






// ------------------------------------------ Optimization Parameters ------------------------------------------- //
const bool HARD_CUTOFF = true;
const bool ENABLE_POTE = true;
const bool IGNORE_HYDROGEN = false;
const int GRIDNODE_QUERY_RANGE = 2;

// If we go larger, a single compound can stretch over 2 nm!
constexpr int MAX_COMPOUND_PARTICLES = 32;
const int MAX_COMPOUNDS = UINT16_MAX-1;			// Arbitrary i think. true max int16_t max - 1. Can also cause trouble when the bondedparticlesLUT static array becomes very large bytewise..

const bool USE_ATOMICS_FOR_BONDS_RESULTS = false;

const int MAX_SAFE_SHIFT = 6;	// Maxmimum manhattan dist that it is safe to shift
