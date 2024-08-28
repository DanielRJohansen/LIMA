#pragma once

#include "Simulation.cuh"
#include "Rasterizer.cuh"
#include "LimaTypes.cuh"
#include "Utilities.h"
#include "MoleculeHull.cuh"

#include <chrono>
#include <string>



class DrawBoxOutlineShader;
class DrawTrianglesShader;
class DrawAtomsShader;
class DrawNormalsShader;

class GLFWwindow;

class FPS {
	std::array<std::chrono::high_resolution_clock::time_point, 20> prevTimepoints;
	int head = 0;
public:
	FPS();
	void NewFrame();
	int GetFps() const;
};


struct FacetTask {
	std::vector<Facet> facets;
	std::optional<Float3> facetsColor;
	bool drawFacetNormals;
	FacetDrawMode facetDrawMode;
};

struct PointsTask {
	std::vector<Float3> points;
	std::optional<Float3> pointsColor;
};




class Display {
public:
	Display(EnvMode);
	~Display();
	void render(const Float3* positions, const std::vector<Compound>& compounds, 
		const BoxParams& boxparams, int64_t step, float temperature, ColoringMethod coloringMethod);	

	void RenderLoop(const MoleculeHullCollection& molCollection, Float3 boxSize, 
		std::optional<std::chrono::milliseconds> duration = std::nullopt);

	// This is meant as a debugging tool
	void RenderLoop(std::vector<FacetTask>& facetTasks, std::vector<PointsTask>& pointsTasks, 
		Float3 boxSize,	std::optional<std::chrono::milliseconds> duration=std::nullopt);

	bool checkWindowStatus();		// Returns false if the windows should close

	int debugValue = 0;

private:
	LimaLogger logger;
	Rasterizer rasterizer;
	bool initGLFW();

	void updateCamera(float pitch, float yaw, float delta_dist=0.f);


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




	std::unique_ptr<DrawBoxOutlineShader> drawBoxOutlineShader = nullptr;
	std::unique_ptr<DrawTrianglesShader> drawTrianglesShader = nullptr;
	std::unique_ptr<DrawAtomsShader> drawAtomsShader = nullptr;
	std::unique_ptr<DrawNormalsShader> drawNormalsShader = nullptr;

	cudaGraphicsResource* renderAtomsBufferCudaResource = nullptr;

	float camera_pitch = 0.f;
	float camera_yaw = 0.f;
	float camera_distance = -2.f;
	Float3 camera_normal{ 0.f,1.f,0.f };

	const std::string window_title = "LIMA - Molecular Dynamics Engine";

	GLFWwindow* window = nullptr;

	const float PI = 3.1415f;

	const int screenHeight = 1400;
	const int screenWidth = 1400;

	const int screensize[2] = {3840, 2160};
};