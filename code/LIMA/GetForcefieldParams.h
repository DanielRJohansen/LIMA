#pragma once

#include "Programs.h"


void GetForcefieldParams() {
	GroFile grofile{ std::filesystem::current_path() / "molecule/conf.gro" };
	TopologyFile topfile{ std::filesystem::current_path() / "molecule/topol.top" };

	Programs::GetForcefieldParams(grofile, topfile, std::filesystem::current_path());
}