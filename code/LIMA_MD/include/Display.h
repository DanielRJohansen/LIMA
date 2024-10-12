#pragma once

#include "Simulation.cuh"
#include "LimaTypes.cuh"
#include "Utilities.h"
#include "MoleculeHull.cuh"

#include <chrono>
#include <string>
#include <thread>
#include <variant>
#include <mutex>
#include <condition_variable>

class DrawBoxOutlineShader;
class DrawTrianglesShader;
class DrawAtomsShader;
class DrawNormalsShader;

class Camera;

class GLFWwindow;

class FPS {
	std::array<std::chrono::high_resolution_clock::time_point, 20> prevTimepoints;
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


//struct FacetTask {
//	std::vector<Facet> facets;
//	std::optional<Float3> facetsColor;
//	bool drawFacetNormals;
//	FacetDrawMode facetDrawMode;
//};
//
//struct PointsTask {
//	std::vector<Float3> points;
//	std::optional<Float3> pointsColor;
//};

namespace Rendering {
	struct SimulationTask {
		const Float3* positions;
		const std::vector<Compound> compounds;
		const BoxParams boxparams;

		std::string siminfo; // Will be output in the window header
		/*int64_t step;
		float temperature;*/
		ColoringMethod coloringMethod;
	};

	struct MoleculehullTask {
		const MoleculeHullCollection& molCollection;
		Float3 boxSize{};
	};

	struct GrofileTask {
		const GroFile& grofile;
		ColoringMethod coloringMethod = Atomname;
	};

	using Task = std::variant<void*, std::unique_ptr<SimulationTask>, std::unique_ptr<MoleculehullTask>, std::unique_ptr<GrofileTask>>;
}




class Display {
public:
	// Functions called by main thread only
	Display(EnvMode);
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

private:
	LimaLogger logger;

	// The renderThread will be spawned during construction, and run this indefinitely
	void Mainloop();

	void Setup();

	bool initGLFW();

	void _RenderAtomsFromCudaresource(Float3 boxSize, int totalParticles);
	void _Render(const MoleculeHullCollection& molCollection, Float3 boxSize);

	void PrepareTask(Rendering::Task& task);

	void PrepareNewRenderTask(const Rendering::SimulationTask&);
	void PrepareNewRenderTask(const Rendering::MoleculehullTask&);
	void PrepareNewRenderTask(const Rendering::GrofileTask&);


	// Interfacing
	bool isDragging = false;
	double lastX = 0.0, lastY = 0.0;
	void OnMouseMove(double xpos, double ypos);
	void OnMouseButton(int button, int action, int mods);
	void OnMouseScroll(double xoffset, double yoffset);
	
	bool pause = false;
	bool renderAtoms = true;
	bool renderFacets = true;
	bool renderFacetsNormals = false;
	FPS fps{};

	Rendering::Task incomingRenderTask = nullptr;
	std::mutex incomingRenderTaskMutex;


	std::unique_ptr<DrawBoxOutlineShader> drawBoxOutlineShader;
	std::unique_ptr<DrawTrianglesShader> drawTrianglesShader;
	std::unique_ptr<DrawAtomsShader> drawAtomsShader;
	std::unique_ptr<DrawNormalsShader> drawNormalsShader;

	cudaGraphicsResource* renderAtomsBufferCudaResource = nullptr;

	std::vector<RenderAtom> renderAtomsTemp;


	std::jthread renderThread;
	std::mutex mutex_;
	std::condition_variable cv_;
	bool setupCompleted = false;

	Camera camera;

	const std::string window_title = "LIMA - Molecular Dynamics Engine";

	GLFWwindow* window = nullptr;

	const float PI = 3.1415f;

	const int screenHeight = 1400;
	const int screenWidth = 1400;

	const int screensize[2] = {3840, 2160};

	std::atomic_bool kill = false;
	std::atomic_bool displaySelfTerminated = false;
};
