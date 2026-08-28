#include "Programs.h"

void SelfTest() {
	const std::filesystem::path workDir = std::filesystem::current_path() / "selftest";

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
		lipidselection.emplace_back(Lipids::Select{ lipidname, workDir, 100. / static_cast<double>(targets.size()) });
	}

	GroFile gro;
	gro.box_size = Float3{ 10.f };
	gro.title = "Membrane";
	TopologyFile top;
	top.SetSystem("Membrane");
	SimulationBuilder::CreateMembrane(gro, top, lipidselection, 5.f);
	Programs::EnergyMinimize(gro, top, false, workDir, Full, true, 5000.f);

	printf("Selftest successful"); // Otherwise we'd have thrown by now
}
