#pragma once

#include "Simulation.cuh"
#include "Rasterizer.cuh"
#include "LimaTypes.cuh"
#include "Utilities.h"

#include <chrono>
#include <string>

#include <GL/glew.h>
#include <glfw3.h>

#ifndef __linux__


const float deg2rad = 2.f * PI / 360.f;
const float rad2deg = 1.f / deg2rad;

class Display {
public:
	Display(EnvMode);
	~Display();
	void render(const Float3* positions, const std::vector<Compound>& compounds, 
		const BoxParams& boxparams, int64_t step, float temperature, ColoringMethod coloringMethod);	

	bool checkWindowStatus();		// Returns false if the windows should close

private:
	LimaLogger logger;
	Rasterizer rasterizer;
	bool initGLFW();

	void updateCamera(float pitch, float yaw, float delta_dist=0.f);



	// OpenGL functions
	void initializePipeline(size_t numAtoms);
	void DrawAtoms(size_t numAtoms);
	void DrawBoxOutline();
	void TerminateGLEW();
	std::optional<GLuint> drawBoxShaderProgram = std::nullopt;
	std::optional<GLuint> drawAtomsShaderProgram = std::nullopt;
	std::optional<GLuint> renderAtomsBuffer = std::nullopt;
	cudaGraphicsResource* renderAtomsBufferCudaResource;
	bool pipelineInitialized = false;


	float camera_pitch = 0.f;
	float camera_yaw = 0.f;
	float camera_distance = -2.f;
	Float3 camera_normal{ 0.f,1.f,0.f };

	const std::string window_title = "LIMA - Molecular Dynamics Engine";

	GLFWwindow* window = nullptr;

	const int triangleAmount = 10; //# of triangles used to draw circle
	const float PI = 3.1415f;
	const float twicePi = 2.0f * PI;

	const int screenHeight = 1400;
	const int screenWidth = 1400;

	const int screensize[2] = {3840, 2160};




};

#else

class Display {
public:
		Display(EnvMode){};
		void render(const Float3* positions, const std::vector<Compound>& compounds, 
		const BoxParams& boxparams, int64_t step, float temperature, ColoringMethod coloringMethod) {}
		bool checkWindowStatus() {return true;}
		void terminate() {}

private:
	int height=-1;
};


#endif


