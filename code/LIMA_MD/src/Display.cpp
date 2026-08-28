#include <GL/glew.h>
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "backends/imgui_impl_opengl3.h"


#include "Display.h"
#include "Shaders.h"    
#include "TimeIt.h"
#include "MDFiles.h"
#include "SSBO.h"



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

void Display::SetupCallbacks() {

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
                case GLFW_KEY_P: {
                    std::lock_guard<std::mutex> lock2(display->liveEditCommandsQueueMutex);
                    display->liveEditCommandsQueue.push_back(LiveEdit::TogglePause{});
                    break;
                }
                case GLFW_KEY_S: {
                    std::lock_guard<std::mutex> lock2(display->liveEditCommandsQueueMutex);
                    display->liveEditCommandsQueue.push_back(LiveEdit::StepOnce{});
					break;
                }
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
}


void Display::Setup() {
    SetThreadName("RenderThread");
    int success = initGLFW();
    SetWindowIcon(window, (FileUtils::GetLimaDir() / "resources"/"logo" / "Lima_Symbol_64x64.png").string().c_str());

   
    // Initialize GLEW
    glewExperimental = GL_TRUE; // Ensure GLEW uses modern techniques for managing OpenGL functionality
    if (glewInit() != GLEW_OK) {
        std::cerr << "Failed to initialize GLEW" << std::endl;
    }

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glfwSetWindowUserPointer(window, this);

    SetupCallbacks();

    overlay = std::make_unique<Overlay>(window, FileUtils::GetLimaDir());

    renderAtomsBuffer = std::make_unique<SSBO>();

    {
        std::lock_guard<std::mutex> lock(mutex_);
        setupCompleted = true; // Set the flag to true after setup is complete
        cv_.notify_one();
    }
}

Display::Display() :
    camera(Float3{ 2.f })
{
        // todo: Display needs to know if its in liveedit, so it knows whether to render gizmo and console.
        // It should also tell Overlay if there is any solvent present, if not dont have the button there..
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



void Display::PrepareTask(Task& task, bool ignorePosition) {
    std::visit([&](auto&& taskPtr) {
        using T = std::decay_t<decltype(taskPtr)>;
        if constexpr (std::is_same_v<T, std::unique_ptr<SimulationTask>>) {
            PrepareNewRenderTask(*taskPtr, ignorePosition);
        }
        else if constexpr (std::is_same_v<T, std::unique_ptr<MoleculehullTask>>) {
            PrepareNewRenderTask(*taskPtr);
        }
        else if constexpr(std::is_same_v<T, std::unique_ptr<GrofileTask>>) {
			PrepareNewRenderTask(*taskPtr);
		}
		else {
			throw std::runtime_error("Unknown task type");
		}
        }, task);
}

void Display::Mainloop() {
    Rendering::Task currentRenderTask = Rendering::NoTask{};
    
    TimeIt frameTime{};

    while (!kill) {
        // Update camera, check if window is closed
        glfwPollEvents();
        if (glfwWindowShouldClose(window)) {
            break;
            printf("Window closed");
        }
        
        bool shouldRecolorAtoms = false;
        ConsumeInputs(shouldRecolorAtoms);

        // Check for new task
        bool newTask = false;
        bool updatedPositions = false;
        {
			std::lock_guard<std::mutex> lock(incomingRenderTaskMutex);     
            if (!incomingRenderTasks.empty()) {
                Rendering::Task incomingRenderTask = std::move(incomingRenderTasks.front());
                incomingRenderTasks.pop_front();

                if (std::holds_alternative<std::unique_ptr<SimulationTaskUpdate>>(incomingRenderTask)) {
                    if (std::holds_alternative<std::unique_ptr<SimulationTask>>(currentRenderTask)) {
                        if (std::get<std::unique_ptr<SimulationTaskUpdate>>(incomingRenderTask) == nullptr) {
                            int a = 0;
                        }
                        PrepareNewRenderTask(*std::get<std::unique_ptr<SimulationTask>>(currentRenderTask), *std::get<std::unique_ptr<SimulationTaskUpdate>>(incomingRenderTask));
                        //incomingRenderTask = Rendering::NoTask{};
                        updatedPositions = true;
                    }
                    else {
                        // This shouldn't happen
                    }
                }
                else if (!std::holds_alternative<Rendering::NoTask>(incomingRenderTask)) {
                    currentRenderTask = std::move(incomingRenderTask);
                    //incomingRenderTasks = Rendering::NoTask{};
                    newTask = true;
                }
            }
        }
        if (newTask || shouldRecolorAtoms) {
            bool ignorePosition = !newTask;
            PrepareTask(currentRenderTask, ignorePosition);
        }
        
        // Check for new input
        bool newInput = false;
		std::optional<std::set<int>> newSelection;
        {
			std::lock_guard<std::mutex> lock(inputMutex);
            newSelection = std::exchange(newSelectionInput, std::nullopt);
        }
        if (newSelection.has_value()) {
            _UpdateSelection(newSelection.value());
            newInput = true;
		}



        const int msPerFrame = std::floor(1. / 60. * 1000.);
        bool shouldDraw = newTask || updatedPositions || newInput || frameTime.elapsed().count() > msPerFrame;

        if (shouldDraw) {
            _Render(currentRenderTask);

            fps.NewFrame();
            frameTime = TimeIt{};
        }
    }
}

void Display::Render(Rendering::Task task, bool blocking) {
    {
        if (std::holds_alternative<std::unique_ptr<Rendering::SimulationTaskUpdate>>(task)) {
            if (std::get<std::unique_ptr<SimulationTaskUpdate>>(task) == nullptr) {
                int a = 0;
            }
        }
		std::lock_guard<std::mutex> lock(incomingRenderTaskMutex);

        if (incomingRenderTasks.size() < 20) // With too many tasks, drop incoming
            incomingRenderTasks.push_back(std::move(task));
    }

    if (blocking) {
        while (1) {
            if (debugValue || displaySelfTerminated) {
                debugValue = 0;
                return;
            }
        }
    }
}

void Display::UpdateSelection(const std::set<int>& selection) {

	std::lock_guard<std::mutex> lock(inputMutex);
	newSelectionInput = selection;
}








bool Display::initGLFW() {
    // Initialize the library
    if (!glfwInit()) {
        throw std::runtime_error("\nGLFW failed to initialize");
    }

    GLFWmonitor* primaryMonitor = glfwGetPrimaryMonitor();
    if (!primaryMonitor) {
        glfwTerminate();
        return -1;
    }

    const GLFWvidmode* mode = glfwGetVideoMode(primaryMonitor);
    int displayWidth = mode->width;
    int displayHeight = mode->height;

    windowSize = glm::ivec2{ (int)((float)displayHeight * 0.8f), (int)((float)displayHeight * 0.8f) };


    // Create a windowed mode window and its OpenGL context
    glfwWindowHint(GLFW_FOCUS_ON_SHOW, GLFW_FALSE); // Do not focus the window on creation
    window = glfwCreateWindow(windowSize.x, windowSize.y, window_title.c_str(), NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        return 0;
    }
#ifndef __linux__
    glfwSetWindowPos(window, displayWidth - windowSize.x - 50, 50);
#endif

    // Make the window's context current
    glfwMakeContextCurrent(window);
    return 1;
}

Float3 Convert(const glm::vec3& v) {
    return Float3{ v.x, v.y, v.z };
}

std::optional<LiveEdit::Command> Display::GetLiveEditCommand() {
    if (activeGizmo && (activeGizmo->pullForce || activeGizmo->rotateForce)) {        
        return LiveEdit::MoveMolecule(Convert(activeGizmo->pullForce.value_or(glm::vec3{})), Convert(activeGizmo->rotateForce.value_or(glm::vec3{})));
    }
    if (bool stopMove = stopMovingLiveeditCmd.exchange(false)) {
		return LiveEdit::MoveMolecule{};
    }
    {
        std::lock_guard<std::mutex> lock2(liveEditCommandsQueueMutex);
        if (!liveEditCommandsQueue.empty()) {
            LiveEdit::Command cmd = liveEditCommandsQueue.front();
            liveEditCommandsQueue.pop_front();
            return cmd;
        }
    }
    return std::nullopt;
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
    const auto avgFrameTime = elapsed / (prevTimepoints.size()-1);
    return static_cast<int>(1e9 / avgFrameTime.count());
}

void Display::TestDisplay() {
	Display display{};
	
	const auto position = std::make_unique<Float3>(0.5f, 0.5f, 0.5f);
	BoxParams params;
    params.boxSize = { 3, 2, 1 };
	params.totalParticles = 1;
    std::vector<PersistentCluster> pclusters(1);
    std::vector<PersistentClusterMeta> pcMetas(1);
    pcMetas.front().particleIdsGlobal[0] = 0;
    pcMetas.front().atomLetter[0] = 'l';
	display.Render(std::make_unique<Rendering::SimulationTask>(pclusters, pcMetas, params), true);
	display.Render(std::make_unique<Rendering::SimulationTaskUpdate>(position.get(), nullptr, SimStatus{}), true);
}