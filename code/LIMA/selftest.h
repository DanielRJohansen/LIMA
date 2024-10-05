#include "Environment.h"
#include "Programs.h"

#include <iostream>


void SelfTest() {
	const std::filesystem::path work_dir = std::filesystem::current_path() / "selftest";

	const fs::path slipidsPath = Filehandler::GetLimaDir() / "resources/Slipids";
	std::vector<std::string> targets;
	for (const auto& entry : fs::directory_iterator(slipidsPath)) {
		if (entry.path().extension() == ".gro") {
			std::string base_name = entry.path().stem().string();
			if (fs::exists(slipidsPath / (base_name + ".itp"))) {
				targets.push_back(base_name);
			}
		}
	}

	Lipids::Selection lipidselection;
	for (const auto& lipidname : targets) {
		lipidselection.emplace_back(Lipids::Select{ lipidname, work_dir, 100. / static_cast<double>(targets.size()) });
	}
	auto [gro, top] = Programs::CreateMembrane(work_dir, lipidselection, Float3{ 10.f }, 5.f, Full);

	printf("Selftest successful"); // Otherwise we'd have thrown by now
}