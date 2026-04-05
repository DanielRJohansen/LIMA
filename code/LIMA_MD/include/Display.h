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

class Camera;
class GLFWwindow;

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

	glm::mat4 View();
	glm::mat4 Projection();
	glm::mat4 ViewProjection();

};


struct Arrow {
	glm::vec3 direction = glm::vec3(1.f, 0.f, 0.f);
	std::vector<Vertex> vertices;
	glm::vec4 color;
	Arrow(glm::vec3 direction, glm::vec4 color);
	void Draw(DrawTrianglesShader*, const glm::mat4& MVP,const glm::vec3& position) const;
};
namespace Rendering {
	struct SimulationTask {
		const Float3* positions;
		std::vector<PersistentCluster> pclusters;
		std::vector<PersistentClusterMeta> pcMeta; // TODO: This could just be a ref, since it remains constant?
		const BoxParams boxparams;

		const std::string siminfo; // Will be output in the window header
		ColoringMethod coloringMethod{};
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

	using Task = std::variant<void*, std::unique_ptr<SimulationTask>, std::unique_ptr<MoleculehullTask>, std::unique_ptr<GrofileTask>>;
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
	std::optional<int> activeAxis = std::nullopt;
	bool isDragging = false;

	glm::vec3 dragStartPosition{};
	glm::vec3 dragStartHitPoint{};

	Arrow arrowX{ glm::vec3(1.f, 0.f, 0.f), glm::vec4(1.f, 0.f, 0.f, 1.f) };
	Arrow arrowY{ glm::vec3(0.f, 1.f, 0.f), glm::vec4(0.f, 1.f, 0.f, 1.f) };
	Arrow arrowZ{ glm::vec3(0.f, 0.f, 1.f), glm::vec4(0.f, 0.f, 1.f, 1.f) };

	void Draw(DrawTrianglesShader* shader, const glm::mat4& VP) const;
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

	void PrepareTask(Rendering::Task& task);

	void PrepareNewRenderTask(const Rendering::SimulationTask&);
	void PrepareNewRenderTask(const Rendering::MoleculehullTask&);
	void PrepareNewRenderTask(Rendering::GrofileTask&);


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
	void HandleGizmo(int atomId);

	bool pause = false;
	bool renderAtoms = true;
	bool renderFacets = true;
	bool renderFacetsNormals = false;
	RenderSettings rendersettings;
	FPS fps{};

	std::optional<TranslateGizmo> activeGizmo;

	Rendering::Task incomingRenderTask = nullptr;
	std::mutex incomingRenderTaskMutex;


	std::unique_ptr<DrawBoxOutlineShader> drawBoxOutlineShader;
	std::unique_ptr<DrawFacetsShader> drawFacetsShader;
	std::unique_ptr<DrawAtomsShader<true>> drawAtomsFromCudaShader;
	std::unique_ptr<DrawAtomsShader<false>> drawAtomsFromCpuShader;
	std::unique_ptr<DrawNormalsShader> drawNormalsShader;
	std::unique_ptr<DrawTrianglesShader> drawTrianglesShader;

	cudaGraphicsResource* renderAtomsBufferCudaResource = nullptr;

	std::vector<RenderAtom> renderAtomsTemp;


	std::jthread renderThread;
	std::mutex mutex_;
	std::condition_variable cv_;
	bool setupCompleted = false;

	std::unique_ptr<Overlay> overlay;
	Camera camera;

	const std::string window_title = "LIMA - Molecular Dynamics Engine";

	GLFWwindow* window = nullptr;
	int2 windowSize{};

	const float PI = 3.1415f;

	std::atomic_bool kill = false;
	std::atomic_bool displaySelfTerminated = false;
};


