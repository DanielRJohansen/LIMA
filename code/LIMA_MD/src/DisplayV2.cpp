#include "DisplayV2.h"
#include "Shaders.h"    

#include <GL/glew.h>
#include <GLFW/glfw3.h>





void Display::updateCamera(float pitch, float yaw, float delta_distance) {
    camera_distance += delta_distance;  
    camera_pitch = pitch;
    camera_yaw = yaw;
    camera_normal = Float3{
    sin(yaw) * cos(pitch),
    cos(yaw) * cos(pitch),
    -sin(pitch)
    };
}

Display::Display(EnvMode envmode) :
    logger(LimaLogger::LogMode::compact, envmode, "display") 
{    
    int success = initGLFW();



    // Initialize GLEW
    glewExperimental = GL_TRUE; // Ensure GLEW uses modern techniques for managing OpenGL functionality
    if (glewInit() != GLEW_OK) {
        std::cerr << "Failed to initialize GLEW" << std::endl;
    }

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glfwSetWindowUserPointer(window, this);

    auto keyCallback = [](GLFWwindow* window, int key, int scancode, int action, int mods) {
        if (action == GLFW_PRESS) {
            // Retrieve the Display instance from the window user pointer
            Display* display = static_cast<Display*>(glfwGetWindowUserPointer(window));
            if (display) {
                const float delta = 3.1415 / 8.f;
                switch (key) {
                case GLFW_KEY_UP:
                    display->updateCamera(display->camera_pitch + delta, display->camera_yaw);
                    break;
                case GLFW_KEY_DOWN:
                    display->updateCamera(display->camera_pitch - delta, display->camera_yaw);
                    break;
                case GLFW_KEY_LEFT:
                    display->updateCamera(display->camera_pitch, display->camera_yaw + delta);
                    break;
                case GLFW_KEY_RIGHT:
                    display->updateCamera(display->camera_pitch, display->camera_yaw - delta);
                    break;
                case GLFW_KEY_PAGE_UP:
                    display->updateCamera(display->camera_pitch, display->camera_yaw, 0.5f);
                    break;
                case GLFW_KEY_PAGE_DOWN:
					display->updateCamera(display->camera_pitch, display->camera_yaw, -0.5f);
					break;
                case GLFW_KEY_N:
					display->debugValue = 1;
					break;
                case GLFW_KEY_P:
                    display->pause = !display->pause;
                    break;
                case GLFW_KEY_1:
                    display->renderAtoms = !display->renderAtoms;
                    break;
                case GLFW_KEY_2:
                    display->renderFacets = !display->renderFacets;
					break;
                }                              
            }
        }
    };
    glfwSetKeyCallback(window, keyCallback);


    glfwSetCursorPosCallback(window, [](GLFWwindow* window, double xpos, double ypos) {
        Display* display = static_cast<Display*>(glfwGetWindowUserPointer(window));
        if (display) {
            display->OnMouseMove(xpos, ypos);
        }
        });

    glfwSetMouseButtonCallback(window, [](GLFWwindow* window, int button, int action, int mods) {
        Display* display = static_cast<Display*>(glfwGetWindowUserPointer(window));
        if (display) {
            display->OnMouseButton(button, action, mods);
        }
        });

    glfwSetScrollCallback(window, [](GLFWwindow* window, double xoffset, double yoffset) {
        Display* display = static_cast<Display*>(glfwGetWindowUserPointer(window));
        if (display) {
            display->OnMouseScroll(xoffset, yoffset);
        }
        });


    logger.finishSection("Display initialized");
}

Display::~Display() {
    glfwTerminate();
}




void Display::OnMouseMove(double xpos, double ypos) {
    if (isDragging) {
        float xOffset = static_cast<float>(xpos - lastX);
        float yOffset = static_cast<float>(lastY - ypos); // Reversed since y-coordinates go from bottom to top

        const float sensitivity = 0.001f; // Adjust sensitivity as needed
        xOffset *= sensitivity;
        yOffset *= sensitivity;

        camera_yaw += xOffset;
        camera_pitch -= yOffset;

        updateCamera(camera_pitch, camera_yaw);
    }

    lastX = xpos;
    lastY = ypos;
}

void Display::OnMouseButton(int button, int action, int mods) {
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        if (action == GLFW_PRESS) {
            isDragging = true;
            glfwGetCursorPos(window, &lastX, &lastY);
        }
        else if (action == GLFW_RELEASE) {
            isDragging = false;
        }
    }
}

void Display::OnMouseScroll(double xoffset, double yoffset) {
	camera_distance += yoffset * 0.1f;
}





bool Display::checkWindowStatus() {
    glfwPollEvents();
    if (glfwWindowShouldClose(window)) {
        glfwTerminate();

        return false;
    }

    return true;
}

bool Display::initGLFW() {
    
    logger.print("Initializing display...\n");

    // Initialize the library
    if (!glfwInit()) {
        throw std::runtime_error("\nGLFW failed to initialize");
    }

    // Create a windowed mode window and its OpenGL context
    logger.print("Loading window --->");
    window = glfwCreateWindow(screenWidth, screenHeight, window_title.c_str(), NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        return 0;
    }
    glfwSetWindowPos(window, screensize[0] - screenWidth - 550, 50);
    logger.print("done\n");

    // Make the window's context current
    glfwMakeContextCurrent(window);
    return 1;
}












    // Constructor initializes the start time and the frame counter
FPS::FPS() {
    auto now = std::chrono::high_resolution_clock::now();
    for (auto& timepoint : prevTimepoints) {
		timepoint = now;
	}
}

    // Call this function when a new frame is rendered
void FPS::NewFrame() {
	using namespace std::chrono;
    head = (head + 1) % prevTimepoints.size();
	prevTimepoints[head] = high_resolution_clock::now();
}

// Returns the current FPS value
int FPS::GetFps() const {
	const int back = (head + 1) % prevTimepoints.size();
    const auto elapsed = duration_cast<std::chrono::nanoseconds>(prevTimepoints[head] - prevTimepoints[back]);
    const auto avgFrameTime = elapsed / prevTimepoints.size();
    return static_cast<int>(1e9 / avgFrameTime.count());
}

