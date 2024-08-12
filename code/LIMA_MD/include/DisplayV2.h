#pragma once

#include "Simulation.cuh"
#include "Rasterizer.cuh"
#include "LimaTypes.cuh"
#include "Utilities.h"

#include <chrono>
#include <string>

#include <GL/glew.h>
#include <GLFW/glfw3.h>



const float deg2rad = 2.f * PI / 360.f;
const float rad2deg = 1.f / deg2rad;
class DrawBoxOutlineShader;
class DrawTrianglesShader;

class Display {
public:
	Display(EnvMode);
	~Display();
	void render(const Float3* positions, const std::vector<Compound>& compounds, 
		const BoxParams& boxparams, int64_t step, float temperature, ColoringMethod coloringMethod);	


	

	void Render(const std::vector<MoleculeContainerSmall>& molecules, float boxlenNM);


	bool checkWindowStatus();		// Returns false if the windows should close

private:
	LimaLogger logger;
	Rasterizer rasterizer;
	bool initGLFW();

	void updateCamera(float pitch, float yaw, float delta_dist=0.f);



	std::unique_ptr<DrawBoxOutlineShader> drawBoxOutlineShader = nullptr;
	std::unique_ptr<DrawTrianglesShader> drawTrianglesShader = nullptr;
	
	// OpenGL functions
	void initializePipeline(size_t numAtoms);
	void DrawAtoms(size_t numAtoms);
	//void DrawBoxOutline();
	void DrawMoleculeContainers(const std::vector<MoleculeContainerSmall>& molecules, float boxlenNM);
	void TerminateGLEW();
	std::optional<GLuint> drawBoxShaderProgram = std::nullopt;
	std::optional<GLuint> drawAtomsShaderProgram = std::nullopt;
	//std::optional<GLuint> drawMoleculeContainersProgram = std::nullopt;
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