#pragma once

#include "Simulation.cuh"
#include "LimaTypes.cuh"
#include "Utilities.h"
#include "MoleculeHull.cuh"
#include "filesystem"
#include "LiveEditCommands.h"
#include "RenderCommons.h"

#include <chrono>
#include <string>
#include <thread>
#include <variant>
#include <mutex>
#include <condition_variable>
#include <set>
#include <deque>

class DrawBoxOutlineShader;
class DrawFacetsShader;
class DrawAtomsShader;
class DrawAtomsPrettyShader;
class DrawNormalsShader;
class DrawTrianglesShader;
class DrawBackgroundGradientShader;

class RenderTargetControl;
class Camera;
class GLFWwindow;
class SSBO;

namespace LimaMoleculeGraph {
	class MoleculeGraph;
}

class FPS {
	std::array<std::chrono::high_resolution_clock::time_point, 32> prevTimepoints;
	int head = 0;
public:
	FPS();
	void NewFrame();
	int GetFps() const;
};

class Camera {
	Float3 center;
	float dist = -2.f;
	float yaw = 0;
	float pitch = 0;
	float aspectRatio = 1.f;

public:
	Camera(Float3 boxSize);
	void Update(float deltaYaw, float deltaPitch, float deltaDist);
	void Update(Float3 boxSize);
	void UpdateViewport(glm::ivec2 viewportSize);

	glm::mat4 View() const;
	glm::mat4 Projection() const;
	glm::mat4 ViewProjection() const;

};

//enum AtomColoringMethod { Name, Charge, Force, GlobalParticleId, PcId};

namespace Rendering {
	struct NoTask{};

	// Sent at simulation start
	struct SimulationTask {
		std::vector<PersistentCluster> pclusters;
		std::vector<PersistentClusterMeta> pcMeta; // TODO: This could just be a ref, since it remains constant?
		const BoxParams boxparams;	
		SimStatus simStatus;
	};
	// Sent at each render-step
	struct SimulationTaskUpdate {
		const Float3* const positions = nullptr;
		const float* const forceMagnitudes = nullptr;
		SimStatus simStatus;
	};

	struct MoleculehullTask {
		const MoleculeHullCollection& molCollection;
		Float3 boxSize{};
	};

	struct GrofileTask {
		const GroFile& grofile;
		bool drawSolvent = true;
		int nAtoms;
		std::set<int> highlightedAtoms;
	};

	using Task = std::variant<NoTask, std::unique_ptr<SimulationTask>, std::unique_ptr<SimulationTaskUpdate>, std::unique_ptr<MoleculehullTask>, std::unique_ptr<GrofileTask>>;
}

struct RenderSettings {
	bool showSolvents = true;
	ColoringMethod coloringMethod{};
};

class Overlay {
public:	
	struct SubmittedCmd { std::string cmd{}; };
	using Command = std::variant<SubmittedCmd, ColoringMethod>;
private:
	bool didDrawThisFrame = false;

	void HandleConsole();
	void HandleContextMenu(RenderSettings& renderSettings, std::optional<glm::dvec2> rightClickedPos);

public:
	std::deque<Command> submittedCommands;
	bool enableConsole = false;	


	Overlay(GLFWwindow*, const std::filesystem::path& limadir);
	~Overlay();

	void Draw(RenderSettings&, const SimStatus&, int fps, std::optional<glm::dvec2> rightClickedPos);
	void Render();
};

struct Arrow {
	glm::vec3 direction = glm::vec3(1.f, 0.f, 0.f);
	std::vector<Vertex> vertices;
	glm::vec4 color;
	int uniqueId;
	Arrow(glm::vec3 direction, glm::vec4 color, int uniqueId);
	void Draw(DrawTrianglesShader*, const glm::mat4& MVP, const glm::vec3& position, float scale = 1.f) const;
};

struct Ring {
	glm::vec3 normal = glm::vec3(0.f, 0.f, 1.f);
	std::vector<Vertex> vertices;
	glm::vec4 color;
	int uniqueId;

	Ring(glm::vec3 normal, glm::vec4 color, int uniqueId);
	void Draw(DrawTrianglesShader* shader, const glm::mat4& VP, const glm::vec3& position, float scale = 1.f) const;
};
struct TransformGizmo {
	glm::vec3 position{};
	int idOfAtomAttachedTo = -1;

	std::optional<int> activeAxis = std::nullopt;
	enum GizmoMode { Translate, Rotate } activeMode = GizmoMode::Translate;

	std::optional<glm::vec3> pullForce;
	std::optional<glm::vec3> rotateForce;

	glm::vec3 dragStartPosition{};
	glm::vec2 dragStartMousePos{};
	glm::vec3 dragStartRotateVector{};

	Arrow arrowX{ glm::vec3(1.f, 0.f, 0.f), glm::vec4(1.f, 0.f, 0.f, 1.f), (int)UniqueRenderElementIds::gizmoArrowX };
	Arrow arrowY{ glm::vec3(0.f, 1.f, 0.f), glm::vec4(0.f, 1.f, 0.f, 1.f), (int)UniqueRenderElementIds::gizmoArrowY };
	Arrow arrowZ{ glm::vec3(0.f, 0.f, 1.f), glm::vec4(0.f, 0.f, 1.f, 1.f), (int)UniqueRenderElementIds::gizmoArrowZ };

	Ring ringX{ glm::vec3(1.f, 0.f, 0.f), glm::vec4(1.f, 0.25f, 0.25f, 1.f), (int)UniqueRenderElementIds::gizmoRotateX };
	Ring ringY{ glm::vec3(0.f, 1.f, 0.f), glm::vec4(0.25f, 1.f, 0.25f, 1.f), (int)UniqueRenderElementIds::gizmoRotateY };
	Ring ringZ{ glm::vec3(0.f, 0.f, 1.f), glm::vec4(0.25f, 0.25f, 1.f, 1.f), (int)UniqueRenderElementIds::gizmoRotateZ };

	void Draw(DrawTrianglesShader* shader, const glm::mat4& VP) const;
	void SetActiveAxis(int selectedObjectId);
	void BeginDragging(glm::vec2 mousePos, const Camera& camera, glm::vec2 windowSize);
	void UpdateDraggingForce(glm::vec2 mousePos, const Camera& camera, glm::vec2 windowSize);
};


class Display {
public:
	// Functions called by main thread only
	Display();
	~Display();
	void WaitForDisplayReady();

	/// <summary>
	/// Queue up a new task for the Display to render. The function waits for a mutex, transfers the data and then leaves
	/// </summary>
	/// <param name=""></param>
	/// <param name="blocking"> If true, the calling thread will be blocked untill debugvalue is set</param>
	void Render(Rendering::Task, bool blocking=false);
	bool DisplaySelfTerminated() { return displaySelfTerminated; }

	void UpdateSelection(const std::set<int>& particleIds);

	volatile int debugValue = 0;

	std::exception_ptr displayThreadException{ nullptr };

	std::atomic_bool allowUserInputs = false;

	static void TestDisplay();
	static void RenderGrofile(const GroFile& grofile, bool drawSolvent=true) {
		Display d;
		d.Render(std::make_unique<Rendering::GrofileTask>(grofile, drawSolvent), true);
	}

	std::optional<LiveEdit::Command> GetLiveEditCommand();
private:
	// The renderThread will be spawned during construction, and run this indefinitely
	void Mainloop();

	void Setup();
	void SetupCallbacks();
	bool ApplyPendingFramebufferResize();

	bool initGLFW();

	void _RenderAtoms();
	void _Render(const MoleculeHullCollection& molCollection, Float3 boxSize);
	void _Render(const Rendering::Task& currentRenderTask); // Render all the things

	void PrepareTask(Rendering::Task& task, bool ignorePosition);

	void PrepareNewRenderTask(const Rendering::SimulationTask&, bool ignorePosition);
	void PrepareNewRenderTask(Rendering::SimulationTask& currentTask, const Rendering::SimulationTaskUpdate&);
	void PrepareNewRenderTask(const Rendering::MoleculehullTask&);
	void PrepareNewRenderTask(Rendering::GrofileTask&);

	void _UpdateSelection(const std::set<int>& selection);


	// Interfacing
	bool isDragging = false;
	glm::dvec2 mousePosAtBtnDown{};
	std::optional<glm::dvec3> mousePosAtRightBtnDown{};
	std::chrono::time_point<std::chrono::steady_clock> timeAtBtnDown;
	glm::dvec2 mousePos{};
	int lastSelectedAtomId = -1;
	void OnMouseMove(double xpos, double ypos);
	void OnMouseButton(int button, int action, int mods);
	void OnMouseScroll(double xoffset, double yoffset);
	void OnMouseLeft();
	void OnMouseLeftClick();
	void HandleGizmo(int atomId);
	int GetObjectIdAtPixel(glm::ivec2);
	void ConsumeInputs(bool& shouldRecolorAtoms);

	std::mutex liveEditCommandsQueueMutex;
	std::deque<LiveEdit::Command> liveEditCommandsQueue;
	bool renderAtoms = true;
	bool renderFacets = true;
	bool renderFacetsNormals = false;
	RenderSettings rendersettings;
	FPS fps{};

	std::optional<TransformGizmo> activeGizmo;
	std::atomic<bool> stopMovingLiveeditCmd = false;

	// Inputs
	std::mutex incomingRenderTaskMutex;
	std::deque<Rendering::Task> incomingRenderTasks;

	std::mutex inputMutex;
	std::optional<std::set<int>> newSelectionInput;
	//
	

	// Shaders
	std::unique_ptr<DrawBoxOutlineShader> drawBoxOutlineShader;
	std::unique_ptr<DrawFacetsShader> drawFacetsShader;
	std::unique_ptr<DrawAtomsShader> drawAtomsFromCpuShader;
	std::unique_ptr<DrawNormalsShader> drawNormalsShader;
	std::unique_ptr<DrawTrianglesShader> drawTrianglesShader;
	std::unique_ptr<DrawBackgroundGradientShader> drawBackgroundGradientShader;
	std::unique_ptr<DrawAtomsPrettyShader> drawAtomsPrettyShader; 

	// Render Data
	cudaGraphicsResource* renderAtomsBufferCudaResource = nullptr;
	std::vector<RenderAtom> renderAtomsHost;
	std::unique_ptr<SSBO> renderAtomsBuffer;


	std::unique_ptr<RenderTargetControl> renderTargetControl;


	std::jthread renderThread;
	std::mutex mutex_;
	std::condition_variable cv_;
	bool setupCompleted = false;

	std::unique_ptr<Overlay> overlay;
	Camera camera;

	const std::string window_title = "LIMA - Molecular Dynamics Engine";

	GLFWwindow* window = nullptr;
	// Cursor positions use logical window coordinates, while OpenGL resources
	// use framebuffer pixels (which can differ on high-DPI displays).
	glm::ivec2 windowSize{};
	glm::ivec2 framebufferSize{};
	bool framebufferResizePending = false;

	const float PI = 3.1415f;

	std::atomic_bool kill = false;
	std::atomic_bool displaySelfTerminated = false;
};


