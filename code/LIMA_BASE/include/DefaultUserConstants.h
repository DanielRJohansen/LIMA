#pragma once

namespace UserConstants {
	// -------------------------------------------- Physics Parameters ---------------------------------------------- //
	constexpr float CUTOFF_NM = 1.2f;

	// ------------------------------------------------ Box Parameters ---------------------------------------------- //
	//constexpr int _BOX_LEN_PM = 7000;
	//constexpr int _BOX_LEN_PM = 18000;
	//constexpr int _BOX_LEN_PM = 20000;
	constexpr int boxlen = 7;	// I need to somehow move this into a namespace
	// ------------------------------------------- Temperature Parameters ------------------------------------------- //
	constexpr bool APPLY_THERMOSTAT = false;		// Apply scalar based on temp
}