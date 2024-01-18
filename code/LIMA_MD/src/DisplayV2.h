#pragma once

#include "Simulation.cuh"
#include "Rasterizer.cuh"
#include "LimaTypes.cuh"
#include "Utilities.h"

#include <chrono>
#include <string>
#include <glfw3.h>

#ifndef __linux__


class Display {
public:
	Display();
	~Display();
	void render(const Float3* positions, const std::vector<Compound>& compounds, 
		const BoxParams& boxparams, int64_t step, float temperature);
	void animate(Trajectory* traj);

	bool checkWindowStatus();		// Returns false if the windows should close
	void terminate();


	//void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
	//	//const float delta = pi / 4.f;

	//	if (action == GLFW_RELEASE) {
	//		/*     switch (key) {
	//			 case GLFW_KEY_UP:
	//				 updateCamera(camera_pitch - delta, camera_yaw);
	//				 break;
	//			 case GLFW_KEY_DOWN:
	//				 updateCamera(camera_pitch + delta, camera_yaw);
	//				 break;
	//			 case GLFW_KEY_LEFT:
	//				 updateCamera(camera_pitch, camera_yaw - delta);
	//				 break;
	//			 case GLFW_KEY_RIGHT:
	//				 updateCamera(camera_pitch, camera_yaw + delta);
	//				 break;
	//			 }*/
	//	}
	//};


private:
	LimaLogger logger;
	Rasterizer rasterizer;
	bool initGLFW();

	void updateCamera(float pitch, float yaw);

	void drawFilledCircle(const RenderBall& ball);
	void drawBalls(const std::vector<RenderBall>& balls, int n_balls);



	int xyToIndex(int x, int y, int size_x) {
		return (x + y * size_x) * 4;
	}

	float camera_pitch = 0.f;
	float camera_yaw = 0.f;

	Float3 camera_normal{ 0.f,1.f,0.f };

	const std::string window_title = "LIMA - Molecular Dynamics Engine";

	GLFWwindow* window = nullptr;

	const int triangleAmount = 10; //# of triangles used to draw circle
	const float PI = 3.1415f;
	const GLfloat twicePi = 2.0f * PI;

	const int display_height = 1400;
	const int display_width = 1400;

	const int screensize[2] = {3840, 2160};
};

#else

class Display {
public:
		Display(){};
		void render(const Float3* positions, const std::vector<Compound>& compounds,
		const BoxParams& boxparams, int64_t step, float temperature) {}
		bool checkWindowStatus() {return true;}
		void terminate() {}

private:
	int height=-1;
};


#endif


