#pragma once

#include "../../LIMA_MD/display/src/DisplayInternal.h"
#include "../../LIMA_MD/display/src/Shaders.h"
#include "NewCartoonRenderer.h"
#include "LipidRenderer.h"
#include "BoxBuilder.cuh"
#include "RenderDataPipe.h"
#include "imgui.h"
#include <GLFW/glfw3.h>
#include <numeric>
#include <stdexcept>

class DisplayTests {
	static void Require(bool condition, const char* message) {
		if (!condition) {
			std::cerr << message << std::endl;
			throw std::runtime_error(message);
		}
	}

	static void PrepareChanges(Display& display) {
		for (auto& [id, context] : display.renderContexts) {
			if (context.shouldRecolorAtoms) {
				display.PrepareTask(context, context.currentRenderTask, true);
				context.shouldRecolorAtoms = false;
			}
		}
	}

	static glm::ivec2 AtomPixel(Display& display, SimulationId id, int atomId) {
		const auto& context = display.renderContexts.at(id);
		const auto& viewport = display.viewports.at(id);
		const auto& position = context.renderAtomsHost.at(atomId).position;
		const auto clip = context.camera->ViewProjection() * glm::vec4(position.x, position.y, position.z, 1.f);
		return glm::ivec2(viewport.origin + glm::dvec2((clip.x / clip.w + 1.) / 2.,
			(1. - clip.y / clip.w) / 2.) * viewport.size);
	}

public:
	static void Run() {
		std::cout << "Checking tile layouts" << std::endl;
		// Odd sizes and fractional DPI must partition the framebuffer without gaps or overlap.
		for (int count = 1; count <= 32; ++count) {
			std::vector<RenderViewport> tiles;
			for (int i = 0; i < count; ++i) {
				const auto tile = RenderViewport::Tile(i, count, { 803, 601 }, { 1205, 902 }, 43.5);
				Require(tile.pixelSize.x > 0 && tile.pixelSize.y > 0, "Empty tile");
				Require(tile.Contains(tile.origin), "Tile must include its starting edge");
				Require(!tile.Contains(tile.origin + tile.size), "Tile must exclude its ending edge");
				for (const auto& previous : tiles) {
					const auto overlap = glm::min(previous.pixelOrigin + previous.pixelSize, tile.pixelOrigin + tile.pixelSize)
						- glm::max(previous.pixelOrigin, tile.pixelOrigin);
					Require(overlap.x <= 0 || overlap.y <= 0, "Tiles overlap in framebuffer space");
				}
				tiles.push_back(tile);
			}
			for (const auto& tile : tiles)
				Require(std::abs(tile.pixelSize.x - tiles.front().pixelSize.x) <= 1 && std::abs(tile.pixelSize.y - tiles.front().pixelSize.y) <= 1,
					"Tiles differ in size");
			if (count == 1 || count == 4 || count == 8 || count == 16 || count == 32) {
				int area = 0;
				for (const auto& tile : tiles) area += tile.pixelSize.x * tile.pixelSize.y;
				Require(area == 1205 * (902 - tiles.front().pixelOrigin.y), "Grid has framebuffer gaps");
			}
		}
		Require(RenderViewport::Tile(0, 1, {}, {}, 0.).pixelSize == glm::ivec2{}, "Minimized layout is not empty");

		std::cout << "Checking camera fit for narrow tiles" << std::endl;
		Camera tileCamera{ Float3{ 30.f, 5.f, 5.f } };
		tileCamera.UpdateViewport({ 200, 600 });
		for (float x : { 0.f, 30.f }) {
			for (float y : { 0.f, 5.f }) {
				for (float z : { 0.f, 5.f }) {
					const glm::vec4 clip = tileCamera.ViewProjection() * glm::vec4(x, y, z, 1.f);
					Require(std::abs(clip.x) <= clip.w && std::abs(clip.y) <= clip.w,
						"Default camera clipped a box corner in a narrow tile");
				}
			}
		}

		std::cout << "Checking lipid interpretation" << std::endl;
		const fs::path slipids = FileUtils::GetLimaDir() / "resources" / "Slipids";
		struct ExpectedLipid { const char* file; LipidHeadgroup headgroup; std::size_t tails; const char* head; };
		// SITS names its atoms C1..C48 rather than by chain, so it checks the interpretation does not rely on names
		for (const auto& expected : std::initializer_list<ExpectedLipid>{
			{ "POPC", LipidHeadgroup::Choline, 2, "P" }, { "POPE", LipidHeadgroup::Ethanolamine, 2, "P" },
			{ "POPS", LipidHeadgroup::Serine, 2, "P" }, { "DOPG", LipidHeadgroup::Anionic, 2, "P" },
			{ "SM16", LipidHeadgroup::Choline, 2, "P" }, { "SITS", LipidHeadgroup::Choline, 2, "P" },
			{ "cholesterol", LipidHeadgroup::Sterol, 1, "O3" } }) {
			const GroFile grofile{ slipids / (std::string{ expected.file } + ".gro") };
			const LipidModels lipids = BoxBuilder::InterpretLipids(grofile);
			Require(lipids.size() == 1, "Lipid was not recognised as exactly one lipid");
			Require(lipids[0].headgroup == expected.headgroup, "Lipid headgroup misclassified");
			Require(lipids[0].tails.size() == expected.tails, "Lipid has the wrong number of tails");
			Require(grofile.atoms[lipids[0].headParticleId].atomName == expected.head, "Lipid head is the wrong atom");
			for (const auto& tail : lipids[0].tails)
				Require(tail.size() >= 8, "Lipid tail is too short to be a chain");
		}
		Require(BoxBuilder::InterpretLipids(GroFile{ FileUtils::GetLimaDir() / "tests" / "clitests" / "fixtures" / "metsol.gro" }).empty(),
			"A protein in water was interpreted as lipids");

		std::cout << "Checking selection options" << std::endl;
		GroFile mixedSystem{ slipids / "POPC.gro" };
		const int lipidAtoms = static_cast<int>(mixedSystem.atoms.size());
		mixedSystem.atoms.push_back({ 2, SmallString{ "ALA" }, SmallString{ "CA" }, 0, { 1.f, 1.f, 1.f } });
		mixedSystem.atoms.push_back({ 2, SmallString{ "ALA" }, SmallString{ "N" }, 0, { 1.1f, 1.f, 1.f } });
		mixedSystem.atoms.push_back({ 3, SmallString{ "SOL" }, SmallString{ "OW" }, 0, { 2.f, 1.f, 1.f } });
		mixedSystem.atoms.push_back({ 3, SmallString{ "SOL" }, SmallString{ "HW1" }, 0, { 2.1f, 1.f, 1.f } });
		mixedSystem.atoms.push_back({ 4, SmallString{ "NA" }, SmallString{ "NA" }, 0, { 3.f, 1.f, 1.f } });
		{
			Rendering::AtomRenderTask task{ mixedSystem };
			std::vector<int> lipidIds(lipidAtoms);
			std::iota(lipidIds.begin(), lipidIds.end(), 0);
			task.molecules = { { "POPC", 1, 1, lipidIds }, { "Protein", 1, 1, { lipidAtoms, lipidAtoms + 1 } },
				{ "SOL", 1, 1, { lipidAtoms + 2, lipidAtoms + 3 } }, { "NA", 1, 1, { lipidAtoms + 4 } } };
			const auto Find = [](const std::vector<Rendering::SelectionOption>& options, const std::string& label) {
				const auto option = std::ranges::find(options, label, &Rendering::SelectionOption::label);
				return option == options.end() ? nullptr : &*option;
			};

			const int carbon = task.lipidModels.at(0).tails[0].back();
			const auto options = Rendering::SelectionOptions(task, carbon, false);
			Require(Find(options, "This atom (#" + std::to_string(carbon) + ")")->selection.atomIds == std::vector<int>{ carbon }, "This atom selects other atoms");
			for (int id : Find(options, "All C atoms")->selection.atomIds)
				Require(task.atoms[id].atomLetter == 'C', "All C atoms selected a non-carbon");
			Require(Find(options, "This molecule")->selection.atomIds == lipidIds, "This molecule is not the lipid");
			Require(Find(options, "All POPC") != nullptr, "No option for all molecules of the type");
			Require(Find(options, "All lipids")->selection.atomIds == lipidIds, "All lipids is not the lipid");
			Require(Find(options, "All protein")->selection.atomIds == std::vector<int>{ lipidAtoms, lipidAtoms + 1 }, "All protein is wrong");
			Require(Find(options, "All water")->selection.atomIds.size() == 2, "All water is wrong");
			Require(Find(options, "All ions")->selection.atomIds == std::vector<int>{ lipidAtoms + 4 }, "All ions is wrong");
			Require(Find(options, "All other") == nullptr, "A category with no atoms was offered");
			Require(Find(options, "All atoms")->selection.atomIds.size() == task.atoms.size(), "All atoms is wrong");
			Require(Find(options, "Clear selection") == nullptr, "Clear offered without a selection");

			const auto emptySpace = Rendering::SelectionOptions(task, -1, true);
			Require(emptySpace.front().label == "All protein", "Empty space offers per-atom options");
			Require(Find(emptySpace, "Clear selection") != nullptr, "No way to clear the selection");
			const auto ion = Rendering::SelectionOptions(task, lipidAtoms + 4, false);
			Require(Find(ion, "All I atoms") == nullptr, "Ions offered as an element");

			// A membrane-only system: the category is everything, and must still be offered
			const Rendering::AtomRenderTask lipidsOnly{ GroFile{ slipids / "POPC.gro" } };
			const auto lipidOptions = Rendering::SelectionOptions(lipidsOnly, carbon, false);
			Require(Find(lipidOptions, "All lipids") != nullptr, "All lipids missing when the system is only lipids");
			// Coordinate-only tasks get a molecule per residue, Slipids names the POPC residue POP
			Require(Find(lipidOptions, "All POP") != nullptr && Find(lipidOptions, "This molecule") != nullptr,
				"Coordinate-only task has no molecule options");
		}

		std::cout << "Starting OpenGL fixture" << std::endl;
		{
			// Keep the window and OpenGL context on this thread for deterministic assertions.
			Display display{ false };
			display.Setup();
			glfwHideWindow(display.window);
			std::cout << "Preparing nine render contexts" << std::endl;
			Require(glfwGetCurrentContext() == display.window, "OpenGL context is not current");
			GroFile molecule;
			molecule.box_size = { 5.f, 5.f, 5.f };
			molecule.atoms = {
				{ 1, SmallString{ "PRO" }, SmallString{ "C" }, 1, { 2.f, 2.5f, 2.5f } },
				{ 1, SmallString{ "PRO" }, SmallString{ "C" }, 2, { 3.f, 2.5f, 2.5f } },
				{ 2, SmallString{ "SOL" }, SmallString{ "O" }, 3, { 2.5f, 2.5f, 3.5f } }
			};
			for (int id = 0; id < 9; ++id) {
				auto task = std::make_unique<Rendering::AtomRenderTask>(molecule, true);
				task->simStatus.step = 1200 + id;
				task->simStatus.temperature = 300.f + id;
				task->simStatus.maxForce = 1.23e3f;
				task->simStatus.avgStepTime = .842f;
				task->simStatus.simulationPerformance = 205.23f;
				task->molecules = { { "Protein", 1, 1, { 0, 1 } }, { "Water", 1, 1, { 2 } } };
				task->backboneChains = { { { { 0, SecondaryStructure::Coil }, { 1, SecondaryStructure::Coil } } } };
				auto& context = display.renderContexts[id];
				context.currentRenderTask = std::move(task);
				display.PrepareTask(context, context.currentRenderTask, false);
			}
			display.activeSimulationId = 0;
			display.activeRenderContext = &display.renderContexts.at(0);
			display.tiled = true;
			std::cout << "Drawing tiled frame" << std::endl;
			display.RenderFrame({});
			std::cout << "Checking tile picking" << std::endl;
			Require(display.viewports.size() == 9, "Not all contexts were drawn");
			for (int id = 0; id < 9; ++id) {
				display.activeSimulationId = id;
				display.activeRenderContext = &display.renderContexts.at(id);
				Require(display.GetObjectIdAtPixel(AtomPixel(display, id, 0)) == 0, "Picking returned wrong atom in tile");
				Require(!glIsEnabled(GL_SCISSOR_TEST), "Picking leaked scissor state");
				display.SelectMolecule(0);
				Require(display.activeRenderContext->selectedMolecule->number == 1, "Selection history leaked across tiles");
			}
			for (auto& [id, context] : display.renderContexts) context.selectedMolecule.reset();
			display.ApplyRepresentation(ColoringMethod::Charge);
			PrepareChanges(display);
			for (const auto& [id, context] : display.renderContexts)
				for (const auto& atom : context.renderAtomsHost)
					Require(atom.flags.y == static_cast<unsigned int>(ColoringMethod::Charge), "Global representation missed atoms");
			for (int id : { 1, 7 }) {
				display.activeSimulationId = id;
				display.activeRenderContext = &display.renderContexts.at(id);
				display.SelectMolecule(0);
			}
			display.ApplyRepresentation(ColoringMethod::GradientFromAtomid);
			PrepareChanges(display);
			display.RenderFrame({});
			for (int id = 0; id < 9; ++id) {
				const auto& context = display.renderContexts.at(id);
				const auto expected = id == 1 || id == 7 ? ColoringMethod::GradientFromAtomid : ColoringMethod::Charge;
				Require(context.renderAtomsHost[0].flags.y == static_cast<unsigned int>(expected), "Representation did not follow selections across tiles");
				Require(context.renderAtomsHost[2].flags.y == static_cast<unsigned int>(ColoringMethod::Charge), "Representation changed unselected atoms");
			}
			for (auto& [id, context] : display.renderContexts) context.selectedMolecule.reset();

			display.activeSimulationId = 4;
			display.activeRenderContext = &display.renderContexts.at(4);
			display.SelectMolecule(0);
			display.ApplyRepresentation(ColoringMethod::NewCartoon);
			PrepareChanges(display);
			Require(display.activeRenderContext->newCartoonRenderer->HasGeometry(), "Selected backbone has no geometry");
			Require(display.activeRenderContext->renderAtomsHost[2].flags.y == static_cast<unsigned int>(ColoringMethod::Charge), "Cartoon changed an unselected atom");
			Require(display.renderContexts.at(3).renderAtomsHost[0].flags.y == static_cast<unsigned int>(ColoringMethod::Charge), "Cartoon changed another tile");
			Require(display.GetObjectIdAtPixel(AtomPixel(display, 4, 2)) == 2, "Mixed cartoon/atom picking failed");

			display.overlay->submittedCommands.push_back(Overlay::SolventVisibility{ false });
			display.ConsumeInputs();
			Require(display.GetObjectIdAtPixel(AtomPixel(display, 4, 2)) == -1, "Invisible solvent intercepted picking");
			display.ApplyRepresentation(ColoringMethod::Atomname);
			PrepareChanges(display);
			Require(!display.activeRenderContext->newCartoonRenderer->HasGeometry(), "Cartoon geometry survived switching selected atoms back");

			// Camera commands broadcast to every visible tile.
			display.overlay->submittedCommands.push_back(Overlay::RevolveCamera{ 1 });
			display.ConsumeInputs();
			for (const auto& [id, context] : display.renderContexts)
				Require(context.revolveCamera, "Tiled orbit did not reach every camera");
			display.overlay->submittedCommands.push_back(Overlay::ResetCamera{ 1 });
			display.ConsumeInputs();
			for (const auto& [id, context] : display.renderContexts)
				Require(!context.revolveCamera, "Tiled reset did not stop every orbit");
			display.isDragging = true;
			display.dragSimulationId = 4;
			const auto draggedView = display.renderContexts.at(4).camera->View();
			const auto otherView = display.renderContexts.at(8).camera->View();
			const auto otherTile = display.viewports.at(8);
			display.OnMouseMove(otherTile.origin.x + otherTile.size.x / 2., otherTile.origin.y + otherTile.size.y / 2.);
			Require(display.renderContexts.at(4).camera->View() != draggedView, "Drag did not move its originating camera");
			Require(display.renderContexts.at(8).camera->View() == otherView, "Drag crossing tiles moved another camera");
			display.framebufferResizePending = true;
			display.ApplyPendingFramebufferResize();
			Require(!display.isDragging && !display.dragSimulationId, "Resize did not cancel drag");
			display.overlay->submittedCommands.push_back(Overlay::SetTiled{ false });
			display.ConsumeInputs();
			display.RenderFrame({});
			Require(display.viewports.size() == 1 && display.viewports.contains(4), "Single view did not retain its context");
			Require(display.GetObjectIdAtPixel(AtomPixel(display, 4, 0)) == 0, "Single-view picking regressed");
			Require(glGetError() == GL_NO_ERROR, "OpenGL error during tile tests");

			std::cout << "Checking lipid representation" << std::endl;
			for (auto& [id, context] : display.renderContexts) context.selectedMolecule.reset();
			{
				auto task = std::make_unique<Rendering::AtomRenderTask>(GroFile{ slipids / "POPC.gro" }, true);
				task->representation = ColoringMethod::Lipids;
				Require(task->lipidModels.size() == 1, "Coordinate-only render task found no lipid");
				const int head = task->lipidModels[0].headParticleId;
				const int body = task->lipidModels[0].tails[0].back();
				auto& context = display.renderContexts[50];
				context.currentRenderTask = std::move(task);
				display.PrepareTask(context, context.currentRenderTask, false);
				display.activeSimulationId = 50;
				display.activeRenderContext = &context;
				display.RenderFrame({});
				Require(context.lipidRenderer && context.lipidRenderer->HasGeometry(), "Lipid representation has no tails");
				Require(context.renderAtomsHost[head].position.w == LipidRenderer::headRadius && context.renderAtomsHost[head].color.w > 0.f,
					"Lipid head is not drawn as a bead");
				Require(context.renderAtomsHost[body].color.w == 0.f, "Lipid body atoms are still drawn");
				Require(display.GetObjectIdAtPixel(AtomPixel(display, 50, head)) == head, "Picking the lipid head failed");
				Require(display.GetObjectIdAtPixel(AtomPixel(display, 50, body)) == head, "Picking a lipid tail did not select its lipid");

				display.ApplyRepresentation(ColoringMethod::Atomname);
				PrepareChanges(display);
				Require(!context.lipidRenderer->HasGeometry(), "Lipid tails survived switching back to atoms");
				Require(context.renderAtomsHost[head].position.w != LipidRenderer::headRadius && context.renderAtomsHost[body].color.w > 0.f,
					"Atoms were not restored after the lipid representation");
			}
			Require(glGetError() == GL_NO_ERROR, "OpenGL error during lipid tests");

			std::cout << "Checking menu selections" << std::endl;
			{
				auto& context = display.renderContexts.at(50);
				const auto* task = std::get<std::unique_ptr<Rendering::AtomRenderTask>>(context.currentRenderTask).get();
				const auto options = Rendering::SelectionOptions(*task, 0, false);
				const auto carbons = std::ranges::find(options, std::string{ "All C atoms" }, &Rendering::SelectionOption::label);
				Require(carbons != options.end(), "No element option for a lipid carbon");
				display.overlay->submittedCommands.push_back(Overlay::Select{ 50, carbons->selection });
				display.ConsumeInputs();
				Require(context.selectedMolecule && context.selectedMolecule->name == "All C atoms", "Menu selection was not applied");
				for (std::size_t id = 0; id < context.renderAtomsHost.size(); ++id)
					Require((context.renderAtomsHost[id].flags.x == 1) == (task->atoms[id].atomLetter == 'C'), "Menu selection highlighted the wrong atoms");

				display.ApplyRepresentation(ColoringMethod::Charge);
				PrepareChanges(display);
				for (std::size_t id = 0; id < context.renderAtomsHost.size(); ++id)
					Require((context.renderAtomsHost[id].flags.y == static_cast<unsigned int>(ColoringMethod::Charge)) == (task->atoms[id].atomLetter == 'C'),
						"Representation did not follow the menu selection");


				display.overlay->submittedCommands.push_back(Overlay::Select{ 50, {} });
				display.ConsumeInputs();
				Require(!context.selectedMolecule, "Clear selection did not clear");
				for (const auto& atom : context.renderAtomsHost)
					Require(atom.flags.x == 0, "Clear selection left highlights");
				// Classic style: the old shader must still draw pickable atoms
				display.overlay->submittedCommands.push_back(Overlay::ClassicAtoms{ true });
				display.ConsumeInputs();
				display.RenderFrame({});
				Require(display.GetObjectIdAtPixel(AtomPixel(display, 50, 0)) >= 0, "Classic atoms drew no pickable atom");
				display.overlay->submittedCommands.push_back(Overlay::ClassicAtoms{ false });
				display.ConsumeInputs();
				Require(display.GetObjectIdAtPixel(AtomPixel(display, 50, 0)) >= 0, "Leaving classic atoms drew no pickable atom");
			}
			display.renderContexts.clear();
			display.activeRenderContext = nullptr;
			auto stoppedPipe = std::make_shared<RenderDataPipe>();
			stoppedPipe->Stop();
			display.renderContexts[99].renderDataPipe = stoppedPipe;
			display.activeSimulationId = 99;
			display.activeRenderContext = &display.renderContexts.at(99);
			Require(display.RemoveStoppedRenderContexts(), "Stopped context was not removed");
			Require(!display.renderContexts.contains(99), "Stopped context remained in the display");
		}
		// Exercise the production ownership path too: ImGui/GL must die before the worker exits.
		for (int i = 0; i < 2; ++i) {
			Display threadedDisplay;
			threadedDisplay.WaitForDisplayReady();
		}
		std::cout << "Display tile layout, picking, selection, representation and routing tests passed\n";
	}
};
