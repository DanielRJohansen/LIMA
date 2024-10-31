#include "Display.h"
#include "Shaders.h"    
#include "TimeIt.h"
#include "MDFiles.h"

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <algorithm>

#define STB_IMAGE_IMPLEMENTATION
#include <stb/stb_image.h>
#undef STB_IMAGE_IMPLEMENTATION


#define NOMINMAX
#if defined(_WIN32) || defined(_WIN64)
#include <windows.h>
#elif defined(__linux__) || defined(__APPLE__)
#include <pthread.h>
#endif

using namespace Rendering;

void SetThreadName(const std::string& name) {
#if defined(_WIN32) || defined(_WIN64)
    // Windows 10, version 1607 and later
    auto handle = GetCurrentThread();
    auto wideName = std::wstring(name.begin(), name.end());
    SetThreadDescription(handle, wideName.c_str());
#elif defined(__linux__) || defined(__APPLE__)
    // pthread_setname_np is available on Linux and macOS
    pthread_setname_np(pthread_self(), name.c_str());
#endif
}
void SetWindowIcon(GLFWwindow* window, const char* iconPath) {
    int width, height, channels;
    unsigned char* pixels = stbi_load(iconPath, &width, &height, &channels, 4);
    if (pixels) {
        GLFWimage icon;
        icon.width = width;
        icon.height = height;
        icon.pixels = pixels;
        glfwSetWindowIcon(window, 1, &icon);
        stbi_image_free(pixels);
    }
    else {
        // Handle error if the image fails to load
        fprintf(stderr, "Failed to load icon image\n");
    }
}


void Display::Setup() {
    int success = initGLFW();
    SetWindowIcon(window, (FileUtils::GetLimaDir() / "resources"/"logo" / "Lima_Symbol_64x64.png").string().c_str());

    SetThreadName("RenderThread");

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
                    display->camera.Update(0, delta, 0);
                    break;
                case GLFW_KEY_DOWN:
                    display->camera.Update(0, -delta, 0);
                    break;
                case GLFW_KEY_LEFT:
                    display->camera.Update(delta, 0, 0);
                    break;
                case GLFW_KEY_RIGHT:
                    display->camera.Update(-delta, 0, 0);
                    break;
                case GLFW_KEY_PAGE_UP:
                    display->camera.Update(0, 0, 0.5f);
                    break;
                case GLFW_KEY_PAGE_DOWN:
                    display->camera.Update(0, 0, -0.5f);
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


    {
        std::lock_guard<std::mutex> lock(mutex_);
        setupCompleted = true; // Set the flag to true after setup is complete
        cv_.notify_one();
    }
}

Display::Display() :
    camera(Float3{ 2.f })
{
    renderThread = std::jthread([this] {
        try {
            Setup();
            Mainloop();
        }
        catch(...) {
            displayThreadException = std::current_exception();
        }
        displaySelfTerminated = true;
    }); 
}

Display::~Display() {
    kill = true;
    if (renderThread.joinable())
        renderThread.join();
    glfwTerminate();
}

void Display::WaitForDisplayReady() {
    if (setupCompleted)
		return;

    // Wait for the render thread to complete Setup
    std::unique_lock<std::mutex> lock(mutex_);
    cv_.wait(lock, [this] { return setupCompleted; });
    // The main thread will block here until setupCompleted becomes true
}



void Display::PrepareTask(Task& task) {
    std::visit([&](auto&& taskPtr) {
        using T = std::decay_t<decltype(taskPtr)>;
        if constexpr (std::is_same_v<T, std::unique_ptr<SimulationTask>>) {
            PrepareNewRenderTask(*taskPtr);
        }
        else if constexpr (std::is_same_v<T, std::unique_ptr<MoleculehullTask>>) {
            PrepareNewRenderTask(*taskPtr);
        }
        else if constexpr(std::is_same_v<T, std::unique_ptr<GrofileTask>>) {
			PrepareNewRenderTask(*taskPtr);
		}
        }, task);
}

void Display::Mainloop() {
    Rendering::Task currentRenderTask = nullptr;

    while (!kill) {
        // Update camera, check if window is closed
        glfwPollEvents();
        if (glfwWindowShouldClose(window)) {
            break;
            printf("Window closed");
        }

        // Check if new data
        {
            bool newData = false;
            incomingRenderTaskMutex.lock();
            if (!std::holds_alternative<void*>(incomingRenderTask)) {
                currentRenderTask = std::move(incomingRenderTask);
                incomingRenderTask = nullptr;
                newData = true;
            }
            incomingRenderTaskMutex.unlock();

            if (newData) {
                PrepareTask(currentRenderTask);
            }
        }

        if (!std::holds_alternative<void*>(currentRenderTask)) {            
            std::visit([&](auto& taskPtr) {
                using T = std::decay_t<decltype(taskPtr)>;
                if constexpr (std::is_same_v<T, std::unique_ptr<SimulationTask>>) {
                    _RenderAtomsFromCudaresource(Float3{ taskPtr->boxparams.boxSize }, taskPtr->boxparams.total_particles);
                }
                else if constexpr (std::is_same_v<T, std::unique_ptr<MoleculehullTask>>) {
                    _Render(taskPtr->molCollection, taskPtr->boxSize);
                }
                else if constexpr(std::is_same_v<T, std::unique_ptr<GrofileTask>>) {
                    _RenderAtomsFromCudaresource(taskPtr->grofile.box_size, taskPtr->nAtoms);
				}
                }, currentRenderTask);
        }
    }
}

void Display::Render(Rendering::Task task, bool blocking) {
    incomingRenderTaskMutex.lock();
    incomingRenderTask = std::move(task);
    incomingRenderTaskMutex.unlock();

    if (blocking) {
        while (1) {
            if (debugValue || displaySelfTerminated) {
                debugValue = 0;
                return;
            }
        }
    }
}









void Display::OnMouseMove(double xpos, double ypos) {
    if (isDragging) {
        const float sensitivity = 0.001f; // Adjust sensitivity as needed
        const float xOffset = static_cast<float>(xpos - lastX) * sensitivity;
        const float yOffset = static_cast<float>(lastY - ypos) * sensitivity; // Reversed since y-coordinates go from bottom to top

        camera.Update(xOffset, -yOffset, 0);
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
    camera.Update(0,0,yoffset * 0.1f);
}

bool Display::initGLFW() {
    // Initialize the library
    if (!glfwInit()) {
        throw std::runtime_error("\nGLFW failed to initialize");
    }

    // Create a windowed mode window and its OpenGL context
    window = glfwCreateWindow(screenWidth, screenHeight, window_title.c_str(), NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        return 0;
    }
#ifndef __linux__
    glfwSetWindowPos(window, screensize[0] - screenWidth - 550, 50);
#endif

    // Make the window's context current
    glfwMakeContextCurrent(window);
    return 1;
}


Camera::Camera(Float3 boxSize) : center(boxSize/2.f), dist(-2.0f * boxSize.y) {}
void Camera::Update(float deltaYaw, float deltaPitch, float deltaDist) {
    yaw += deltaYaw;
    pitch += deltaPitch;
    dist += deltaDist + deltaDist * -std::min(dist, 0.f) * 0.5f;
}
void Camera::Update(Float3 boxSize) {
    if (center != boxSize / 2.f)
        *this = Camera(boxSize);
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

void Display::TestDisplay() {
	Display display{};
	
	const auto position = std::make_unique<Float3>(0.5f, 0.5f, 0.5f);
	Compound compound;
	compound.n_particles = 1;
	compound.atomLetters[0] = '_';
	BoxParams params;
	params.boxSize = 2;
	params.total_compound_particles = 1;
	params.total_particles = 1;
	params.total_particles_upperbound = 1;
	display.Render(std::make_unique<Rendering::SimulationTask>(position.get(), std::vector<Compound>{compound}, params, "", Atomname), true);
}

//void Display::RenderGrofile(const GroFile& grofile) {
//    Display d;
//    d.Render(std::make_unique<Rendering::GrofileTask>(grofile), true);
//}