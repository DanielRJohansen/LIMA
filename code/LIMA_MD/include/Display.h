#pragma once

#include "LiveEditCommands.h"
#include "RenderTask.h"

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <deque>
#include <exception>
#include <memory>
#include <mutex>
#include <optional>
#include <set>
#include <string>
#include <thread>
#include <vector>

#include <glm.hpp>

class DrawBoxOutlineShader;
class DrawFacetsShader;
class DrawAtomsShader;
class DrawAtomsPrettyShader;
class DrawNormalsShader;
class DrawTrianglesShader;
class DrawBackgroundGradientShader;
class RenderTargetControl;
class Camera;
class FPS;
class GLFWwindow;
class Overlay;
struct RenderSettings;
class SSBO;
struct TransformGizmo;
namespace NewCartoon { class Renderer; }

class Display {
public:
	Display();
	~Display();
	void WaitForDisplayReady();

	void Render(Rendering::Task, bool blocking = false);
	bool DisplaySelfTerminated() { return displaySelfTerminated; }

	void UpdateSelection(const std::set<int>& particleIds);
	void SetSpinnerVisible(bool visible) { spinnerVisible.store(visible); }

	volatile int debugValue = 0;
	std::exception_ptr displayThreadException{ nullptr };
	std::atomic_bool allowUserInputs = false;

	static void TestDisplay();
	static void RenderGrofile(const GroFile& grofile, bool showSolvents = true) {
		Display display;
		display.Render(std::make_unique<Rendering::AtomRenderTask>(grofile, showSolvents), true);
	}

	std::optional<LiveEdit::Command> GetLiveEditCommand();

private:
	void Mainloop();
	void Setup();
	void SetupCallbacks();
	bool ApplyPendingFramebufferResize();
	bool initGLFW();

	void _RenderAtoms();
	void _Render(const MoleculeHullCollection& molCollection, Float3 boxSize);
	void _Render(const Rendering::Task& currentRenderTask);
	void PrepareTask(Rendering::Task& task, bool ignorePosition);
	void PrepareNewRenderTask(Rendering::AtomRenderTask&, bool ignorePosition);
	void PrepareNewRenderTask(Rendering::AtomRenderTask& currentTask, const Rendering::SimulationTaskUpdate&);
	void PrepareNewRenderTask(const Rendering::MoleculehullTask&);
	void _UpdateSelection(const std::set<int>& selection);

	void OnMouseMove(double xpos, double ypos);
	void OnMouseButton(int button, int action, int mods);
	void OnMouseScroll(double xoffset, double yoffset);
	void OnMouseLeft();
	void OnMouseLeftClick();
	void HandleGizmo(int atomId);
	int GetObjectIdAtPixel(glm::ivec2);
	void ConsumeInputs(bool& shouldRecolorAtoms);

	bool isDragging = false;
	glm::dvec2 mousePosAtBtnDown{};
	std::optional<glm::dvec3> mousePosAtRightBtnDown{};
	std::chrono::time_point<std::chrono::steady_clock> timeAtBtnDown;
	glm::dvec2 mousePos{};
	int lastSelectedAtomId = -1;

	std::mutex liveEditCommandsQueueMutex;
	std::deque<LiveEdit::Command> liveEditCommandsQueue;
	bool renderAtoms = true;
	bool renderFacets = true;
	bool renderFacetsNormals = false;
	std::unique_ptr<RenderSettings> rendersettings;
	std::unique_ptr<FPS> fps;
	std::unique_ptr<TransformGizmo> activeGizmo;
	std::atomic<bool> stopMovingLiveeditCmd = false;

	std::mutex incomingRenderTaskMutex;
	std::deque<Rendering::Task> incomingRenderTasks;
	std::mutex inputMutex;
	std::optional<std::set<int>> newSelectionInput;

	std::unique_ptr<DrawBoxOutlineShader> drawBoxOutlineShader;
	std::unique_ptr<DrawFacetsShader> drawFacetsShader;
	std::unique_ptr<DrawAtomsShader> drawAtomsFromCpuShader;
	std::unique_ptr<DrawNormalsShader> drawNormalsShader;
	std::unique_ptr<DrawTrianglesShader> drawTrianglesShader;
	std::unique_ptr<DrawBackgroundGradientShader> drawBackgroundGradientShader;
	std::unique_ptr<DrawAtomsPrettyShader> drawAtomsPrettyShader;
	std::unique_ptr<NewCartoon::Renderer> newCartoonRenderer;

	cudaGraphicsResource* renderAtomsBufferCudaResource = nullptr;
	std::vector<RenderAtom> renderAtomsHost;
	std::unique_ptr<SSBO> renderAtomsBuffer;
	std::unique_ptr<RenderTargetControl> renderTargetControl;

	std::jthread renderThread;
	std::mutex mutex_;
	std::condition_variable cv_;
	bool setupCompleted = false;

	std::unique_ptr<Overlay> overlay;
	std::unique_ptr<Camera> camera;
	const std::string window_title = "LIMA - Molecular Dynamics Engine";
	GLFWwindow* window = nullptr;
	glm::ivec2 windowSize{};
	glm::ivec2 framebufferSize{};
	bool framebufferResizePending = false;

	const float PI = 3.1415f;
	std::atomic_bool kill = false;
	std::atomic_bool displaySelfTerminated = false;
	std::atomic_bool spinnerVisible = false;
};
