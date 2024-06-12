#pragma once

namespace UserConstants {
	// -------------------------------------------- Physics Parameters ---------------------------------------------- //
	constexpr float CUTOFF_NM = 1.2f;

	// ------------------------------------------------ Box Parameters ---------------------------------------------- //
	constexpr int boxlen = 23.f;	// I need to somehow move this into a namespace
	// ------------------------------------------- Temperature Parameters ------------------------------------------- //
	constexpr bool APPLY_THERMOSTAT = false;		// Apply scalar based on temp
}