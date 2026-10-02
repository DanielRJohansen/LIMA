#pragma once

#include "../../LIMA_MD/src/DisplayInternal.h"
#include "../../LIMA_MD/src/Shaders.h"
#include "NewCartoonRenderer.h"
#include "RenderDataPipe.h"
#include "imgui.h"
#include <GLFW/glfw3.h>
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
		for (int count = 1; count <= 9; ++count) {
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
			if (count == 1 || count == 4 || count == 9) {
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
