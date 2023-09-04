#pragma once

#include "Simulation.cuh"
#include "Rasterizer.cuh"
#include "LimaTypes.cuh"
#include "Utilities.h"

#include <chrono>
#include <string>

#ifndef __linux__
#include <glfw3.h>
#define ENABLE_DISPLAY 
#else
#include <GLFW/glfw3.h>
#endif

class Display {
public:
	Display();
	~Display();
	void render(Simulation* simulation);
	void animate(Trajectory* traj);

	bool checkWindowStatus();		// Returns false if the windows should close
	void terminate(); 

private:
	LimaLogger logger;
	Rasterizer rasterizer;
	bool initGLFW();

	void drawFilledCircle(const RenderBall& ball);
	void drawBalls(const std::vector<RenderBall>& balls, int n_balls);

	int xyToIndex(int x, int y, int size_x) {
		return (x + y * size_x) * 4;
	}


	const std::string window_title = "LIMA - Molecular Dynamics Engine";

	GLFWwindow* window = nullptr;

	const int triangleAmount = 10; //# of triangles used to draw circle
	const float PI = 3.1415f;
	const GLfloat twicePi = 2.0f * PI;

	const int display_height = 1400;
	const int display_width = 1400;

	const int screensize[2] = {3840, 2160};
};
