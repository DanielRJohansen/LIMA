#include "Environment.h"
#include "Programs.h"

#include <iostream>


void SelfTest() {
	const std::filesystem::path work_dir = std::filesystem::current_path() / "selftest";

	Lipids::Selection lipidselection;
	const std::array<std::string, 6> lipids = { "POPC", "POPE", "DDPC", "DMPC", "Cholesterol", "DOPC" };
	for (const auto& lipidname : lipids) {
		lipidselection.emplace_back(Lipids::Select{ lipidname, work_dir, 100./6. });	// 10% of each lipid, except 50% POPC
	}
	auto [gro, top] = Programs::CreateMembrane(work_dir, lipidselection, Float3{ 7.f }, 3.5f, Full);
	printf("Selftest successful"); // Otherwise we'd have thrown by now
}