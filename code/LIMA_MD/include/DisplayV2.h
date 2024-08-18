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

class Display {
public:
	Display(EnvMode);
	~Display();
	void render(const Float3* positions, const std::vector<Compound>& compounds, 
		const BoxParams& boxparams, int64_t step, float temperature, ColoringMethod coloringMethod);	

	void Render(const MoleculeHullCollection& molCollection, Float3 boxSize);

	// This is meant as a debugging tool
	void Render(const std::vector<Facet>& facets, std::optional<Float3> facetsColor, bool drawFacetNormals, FacetDrawMode,
		const std::vector<Float3>& points, std::optional<Float3> pointsColor,
		Float3 boxSize, bool flushGLBuffer, bool swapGLBuffers);

	bool checkWindowStatus();		// Returns false if the windows should close

	int debugValue = 0;

private:
	LimaLogger logger;
	Rasterizer rasterizer;
	bool initGLFW();

	void updateCamera(float pitch, float yaw, float delta_dist=0.f);



	std::unique_ptr<DrawBoxOutlineShader> drawBoxOutlineShader = nullptr;
	std::unique_ptr<DrawTrianglesShader> drawTrianglesShader = nullptr;
	std::unique_ptr<DrawAtomsShader> drawAtomsShader = nullptr;
	std::unique_ptr<DrawNormalsShader> drawNormalsShader = nullptr;

	cudaGraphicsResource* renderAtomsBufferCudaResource;


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