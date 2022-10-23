#pragma once

#include "Simulation.cuh"
#ifndef __linux__
#include <glfw3.h>
#define ENABLE_DISPLAY 
#else
#include <GLFW/glfw3.h>
#endif // !__linux__



#include "Rasterizer.cuh"
#include "LimaTypes.cuh"

#include <chrono>

class DisplayV2 {
public:
	DisplayV2();
	DisplayV2(Simulation* simulation);
	void render(Simulation* simulation);
	void animate(Trajectory* traj);

	bool checkWindowStatus();		// Returns false if the windows should close
	void terminate(); 
	//sf::RenderWindow* window;
	//Raytracer* raytracer;

private:
	Rasterizer rasterizer;
	bool initGLFW();

	void drawFilledCircle(GLfloat x, GLfloat y, GLfloat radius, Int3 color);

	void drawBalls(RenderBall* balls, int n_balls);
	void draw(uint8_t* image);
	void sleep(int ms);
	uint8_t* enhance(uint8_t* im, int from_size);	// doubles image resolution

	int xyToIndex(int x, int y, int size_x) {
		return (x + y * size_x) * 4;
	}




	GLFWwindow* window;

	const int triangleAmount = 20; //# of triangles used to draw circle
	const float PI = 3.14f;
	const GLfloat twicePi = 2.0f * PI;
};

