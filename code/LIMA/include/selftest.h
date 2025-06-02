#include "Programs.h"

void SelfTest() {
	const std::filesystem::path work_dir = std::filesystem::current_path() / "selftest";

	const fs::path slipidsPath = FileUtils::GetLimaDir() / "resources/Slipids";
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

	auto [gro, top] = SimulationBuilder::CreateMembrane(lipidselection, Float3{ 10.f }, 5.f);
	Programs::EnergyMinimize(*gro, *top, false, work_dir, Full, true, 5000.f);

	printf("Selftest successful"); // Otherwise we'd have thrown by now
}