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
template <bool>class DrawAtomsShader;
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

public:
	Camera(Float3 boxSize);
	void Update(float deltaYaw, float deltaPitch, float deltaDist);
	void Update(Float3 boxSize);

	glm::mat4 View() const;
	glm::mat4 Projection() const;
	glm::mat4 ViewProjection() const;

};


struct Arrow {
	glm::vec3 direction = glm::vec3(1.f, 0.f, 0.f);
	std::vector<Vertex> vertices;
	glm::vec4 color;
	int uniqueId;
	Arrow(glm::vec3 direction, glm::vec4 color, int uniqueId);
	void Draw(DrawTrianglesShader*, const glm::mat4& MVP,const glm::vec3& position, float scale = 1.f) const;
};
namespace Rendering {
	struct NoTask{};

	// Sent at simulation start
	struct SimulationTask {
		std::vector<PersistentCluster> pclusters;
		std::vector<PersistentClusterMeta> pcMeta; // TODO: This could just be a ref, since it remains constant?
		const BoxParams boxparams;
		ColoringMethod coloringMethod{};		
		SimStatus simStatus;
	};
	// Sent at each render-step
	struct SimulationTaskUpdate {
		const Float3* positions;
		SimStatus simStatus;
	};

	struct MoleculehullTask {
		const MoleculeHullCollection& molCollection;
		Float3 boxSize{};
	};

	struct GrofileTask {
		const GroFile& grofile;
		bool drawSolvent = true;
		ColoringMethod coloringMethod = Atomname;
		int nAtoms;
		std::set<int> highlightedAtoms;
	};

	using Task = std::variant<NoTask, std::unique_ptr<SimulationTask>, std::unique_ptr<SimulationTaskUpdate>, std::unique_ptr<MoleculehullTask>, std::unique_ptr<GrofileTask>>;
}

struct RenderSettings {
	bool showSolvents = true;
};

class Overlay {
	bool didDrawThisFrame = false;

	void HandleConsole();

public:
	std::mutex consoleMutex;
	std::deque<std::string> submittedCommands;// If we ever access from other than renderthread, well need a mutex`

	Overlay(GLFWwindow*, const std::filesystem::path& limadir);
	~Overlay();

	void Draw(RenderSettings&, const SimStatus&, int fps);
	void Render();
};


struct TranslateGizmo {
	glm::vec3 position{};
	int idOfAtomAttachedTo = -1;
	std::optional<int> activeAxis = std::nullopt;
	std::optional<glm::vec3> pullForce;
	//std::optional<int> hoveredAxis = std::nullopt;

	glm::vec3 dragStartPosition{};
	glm::vec2 dragStartMousePos{};

	Arrow arrowX{ glm::vec3(1.f, 0.f, 0.f), glm::vec4(1.f, 0.f, 0.f, 1.f), (int)UniqueRenderElementIds::gizmoArrowX };
	Arrow arrowY{ glm::vec3(0.f, 1.f, 0.f), glm::vec4(0.f, 1.f, 0.f, 1.f), (int)UniqueRenderElementIds::gizmoArrowY };
	Arrow arrowZ{ glm::vec3(0.f, 0.f, 1.f), glm::vec4(0.f, 0.f, 1.f, 1.f), (int)UniqueRenderElementIds::gizmoArrowZ };

	void Draw(DrawTrianglesShader* shader, const glm::mat4& VP) const;
	void SetActiveAxis(int selectedObjectId);
	void BeginDragging(glm::vec2 mousePos, const Camera& camera);
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

	bool initGLFW();

	void _RenderAtoms(Float3 boxSize, int totalParticles, bool fromCuda);
	void _Render(const MoleculeHullCollection& molCollection, Float3 boxSize);
	void _Render(const Rendering::Task& currentRenderTask); // Render all the things

	void PrepareTask(Rendering::Task& task);

	void PrepareNewRenderTask(const Rendering::SimulationTask&);
	void PrepareNewRenderTask(Rendering::SimulationTask& currentTask, const Rendering::SimulationTaskUpdate&);
	void PrepareNewRenderTask(const Rendering::MoleculehullTask&);
	void PrepareNewRenderTask(Rendering::GrofileTask&);

	void _UpdateSelection(const std::set<int>& selection);


	// Interfacing
	bool isDragging = false;
	glm::dvec2 mousePosAtBtnDown{};
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
	void ConsumeInputs();

	std::mutex liveEditCommandsQueueMutex;
	std::deque<LiveEdit::Command> liveEditCommandsQueue;
	bool renderAtoms = true;
	bool renderFacets = true;
	bool renderFacetsNormals = false;
	RenderSettings rendersettings;
	FPS fps{};

	std::optional<TranslateGizmo> activeGizmo;
	std::atomic<bool> stopMovingLiveeditCmd = false;

	// Inputs
	std::mutex incomingRenderTaskMutex;
	Rendering::Task incomingRenderTask = Rendering::NoTask{};

	std::mutex inputMutex;
	std::optional<std::set<int>> newSelectionInput;
	//
	

	// Shaders
	std::unique_ptr<DrawBoxOutlineShader> drawBoxOutlineShader;
	std::unique_ptr<DrawFacetsShader> drawFacetsShader;
	std::unique_ptr<DrawAtomsShader<true>> drawAtomsFromCudaShader;
	std::unique_ptr<DrawAtomsShader<false>> drawAtomsFromCpuShader;
	std::unique_ptr<DrawNormalsShader> drawNormalsShader;
	std::unique_ptr<DrawTrianglesShader> drawTrianglesShader;
	std::unique_ptr<DrawBackgroundGradientShader> drawBackgroundGradientShader;

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
	glm::ivec2 windowSize{};

	const float PI = 3.1415f;

	std::atomic_bool kill = false;
	std::atomic_bool displaySelfTerminated = false;
};


