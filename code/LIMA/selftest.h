#include "Environment.h"
#include "Programs.h"

#include <iostream>


void SelfTest() {
	const std::filesystem::path work_dir = std::filesystem::current_path() / "selftest";
	Environment env{ work_dir.string(), EnvMode::Full, false };
	env.CreateSimulation(20.f);
	LipidsSelection lipids;
	lipids.emplace_back(LipidSelect{ "POPC", 50 });
	lipids.emplace_back(LipidSelect{ "DMPC", 40 });
	lipids.emplace_back(LipidSelect{ "cholesterol", 10 });
	Programs::CreateMembrane(env, lipids, true, 3.5f, true);
}