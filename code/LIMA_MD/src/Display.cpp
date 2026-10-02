#include <GL/glew.h>
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "backends/imgui_impl_opengl3.h"


#include "Display.h"
#include "DisplayInternal.h"
#include "Shaders.h"
#include "NewCartoonRenderer.h"
#include "TimeIt.h"
#include "MDFiles.h"
#include "SSBO.h"
#include "RenderDataPipe.h"



#include <GLFW/glfw3.h>
#include <algorithm>
#include <format>



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

RenderContext::RenderContext()
	: renderSettings(std::make_unique<RenderSettings>())
	, camera(std::make_unique<Camera>(Float3{ 2.f })) {}

RenderContext::~RenderContext() = default;
RenderContext::RenderContext(RenderContext&&) noexcept = default;
RenderContext& RenderContext::operator=(RenderContext&&) noexcept = default;

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
	glfwSetWindowSizeCallback(window, [](GLFWwindow* window, int width, int height) {
		Display* display = static_cast<Display*>(glfwGetWindowUserPointer(window));
		if (display)
			display->windowSize = glm::ivec2{ width, height };
	});

	glfwSetFramebufferSizeCallback(window, [](GLFWwindow* window, int width, int height) {
		Display* display = static_cast<Display*>(glfwGetWindowUserPointer(window));
		if (display) {
			display->framebufferSize = glm::ivec2{ width, height };
			display->framebufferResizePending = true;
		}
	});

    auto keyCallback = [](GLFWwindow* window, int key, int scancode, int action, int mods) {
        if (ImGui::GetCurrentContext() && ImGui::GetIO().WantCaptureKeyboard)
            return;
        if (action == GLFW_PRESS) {
            // Retrieve the Display instance from the window user pointer
            Display* display = static_cast<Display*>(glfwGetWindowUserPointer(window));
            if (display) {
                const bool hasMouseContext = display->TargetMouseContext();
                const float delta = 3.1415 / 8.f;
                switch (key) {
                case GLFW_KEY_UP:
                    if (hasMouseContext && display->activeRenderContext) display->activeRenderContext->camera->Update(0, delta, 0);
                    break;
                case GLFW_KEY_DOWN:
                    if (hasMouseContext && display->activeRenderContext) display->activeRenderContext->camera->Update(0, -delta, 0);
                    break;
                case GLFW_KEY_LEFT:
                    if (hasMouseContext && display->activeRenderContext) display->activeRenderContext->camera->Update(delta, 0, 0);
                    break;
                case GLFW_KEY_RIGHT:
                    if (hasMouseContext && display->activeRenderContext) display->activeRenderContext->camera->Update(-delta, 0, 0);
                    break;
                case GLFW_KEY_PAGE_UP:
                    if (hasMouseContext && display->activeRenderContext) display->activeRenderContext->camera->Update(0, 0, 0.5f);
                    break;
                case GLFW_KEY_PAGE_DOWN:
                    if (hasMouseContext && display->activeRenderContext) display->activeRenderContext->camera->Update(0, 0, -0.5f);
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
                    if (hasMouseContext && display->activeRenderContext) display->activeRenderContext->renderAtoms = !display->activeRenderContext->renderAtoms;
                    break;
                case GLFW_KEY_2:
                    if (hasMouseContext && display->activeRenderContext) display->activeRenderContext->renderFacets = !display->activeRenderContext->renderFacets;
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
	glfwGetWindowSize(window, &windowSize.x, &windowSize.y);
	glfwGetFramebufferSize(window, &framebufferSize.x, &framebufferSize.y);
	framebufferResizePending = true;
	ApplyPendingFramebufferResize();

    overlay = std::make_unique<Overlay>(window, FileUtils::GetLimaDir());

    //renderAtomsBuffer = std::make_unique<SSBO>();

    {
        std::lock_guard<std::mutex> lock(mutex_);
        setupCompleted = true; // Set the flag to true after setup is complete
        cv_.notify_one();
    }
}

Display::Display() : Display(true) {}

Display::Display(bool startRenderThread)
	: fps(std::make_unique<FPS>())
{
    if (!startRenderThread) return;
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
        ReleaseGraphics();
        displaySelfTerminated = true;
    }); 
}

Display::~Display() {
    kill = true;
    if (renderThread.joinable())
        renderThread.join();
    else
        ReleaseGraphics();
    glfwTerminate();
}

void Display::ReleaseGraphics() {
    // Windows destroys a thread's windows on exit. Release ImGui and GPU resources
    // on the window's owning thread while its OpenGL context is still valid.
    if (window) {
        glfwMakeContextCurrent(window);
        activeRenderContext = nullptr;
        renderContexts.clear();
        overlay.reset();
        renderTargetControl.reset();
        drawBoxOutlineShader.reset();
        drawFacetsShader.reset();
        drawAtomsFromCpuShader.reset();
        drawNormalsShader.reset();
        drawTrianglesShader.reset();
        drawBackgroundGradientShader.reset();
        drawAtomsPrettyShader.reset();
        glfwDestroyWindow(window);
        window = nullptr;
    }
}

void Display::WaitForDisplayReady() {
    if (setupCompleted)
		return;

    // Wait for the render thread to complete Setup
    std::unique_lock<std::mutex> lock(mutex_);
    cv_.wait(lock, [this] { return setupCompleted; });
    // The main thread will block here until setupCompleted becomes true
}



void Display::PrepareTask(RenderContext& renderContext, Task& task, bool ignorePosition) {
    std::visit([&](auto&& taskPtr) {
        using T = std::decay_t<decltype(taskPtr)>;
		if constexpr (std::is_same_v<T, std::unique_ptr<AtomRenderTask>>) {
            PrepareNewRenderTask(renderContext, *taskPtr, ignorePosition);
        }
        else if constexpr (std::is_same_v<T, std::unique_ptr<MoleculehullTask>>) {
            PrepareNewRenderTask(renderContext, *taskPtr);
        }
		else {
			throw std::runtime_error("Unknown task type");
		}
        }, task);
}

void Display::Mainloop() {

    TimeIt frameTime{};

    while (!kill) {
        // Update camera, check if window is closed
        glfwPollEvents();
		const bool framebufferWasResized = ApplyPendingFramebufferResize();
        if (glfwWindowShouldClose(window)) {
            break;
            printf("Window closed");
        }
        
        // Check for new user input
        bool newInput = false;
		std::deque<std::tuple<SimulationId, std::set<int>, std::optional<Rendering::MoleculeInfo>>> newSelections;
        {
            std::lock_guard<std::mutex> lock(inputMutex);
            newSelections.swap(newSelectionInputs);
        }
        for (auto& newSelection : newSelections) {
            if (auto context = renderContexts.find(std::get<0>(newSelection)); context != renderContexts.end()) {
				auto& renderContext = context->second;
                _UpdateSelection(renderContext, std::get<1>(newSelection));
				renderContext.selectedMolecule = std::get<2>(newSelection);
                newInput = true;
			}
        }

        // Check for newly submitted work to display
        {
            std::lock_guard<std::mutex> lock(incomingRenderTaskMutex);
			while (!incomingRenderTasksGlobal.empty()) {
				auto [simId, incomingRenderTask, renderDataPipe, label] = std::move(incomingRenderTasksGlobal.front());
                incomingRenderTasksGlobal.pop_front();
                if (std::holds_alternative<Rendering::FreeTask>(incomingRenderTask)) {
					CancelInteraction();
					viewports.clear();
					if (popupSimulationId == simId) popupSimulationId.reset();
					if (activeSimulationId == simId) {
						activeSimulationId.reset();
						activeRenderContext = nullptr;
					}
					renderContexts.erase(simId);
                }
                else {
					const bool isNewContext = !renderContexts.contains(simId);
                    if (isNewContext) {
                        CancelInteraction();
                        viewports.clear();
                    }
                    auto& context = renderContexts[simId];
					if (renderDataPipe)
						context.renderDataPipe = renderDataPipe;
					if (!label.empty())
						context.label = std::move(label);
                    if (context.incomingRenderTasks.size() < 10)
					    renderContexts[simId].incomingRenderTasks.push_back(std::move(incomingRenderTask));
					if (isNewContext && renderContexts.size() > 1 && !allowUserInputs)
						tiled = true;
                }
            }
        }

		RemoveStoppedRenderContexts();
		if (renderContexts.size() <= 1)
			tiled = false;


        if (renderContexts.empty()) {
            activeRenderContext = nullptr;
            viewports.clear();
            if (framebufferSize.x > 0 && framebufferSize.y > 0) {
				overlay->enableConsole = allowUserInputs;
				RenderSettings menuSettings{};
				overlay->BeginFrame(menuSettings, fps->GetFps(), {}, false, 0);
				glDisable(GL_SCISSOR_TEST);
				glViewport(0, 0, framebufferSize.x, framebufferSize.y);
                glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
				overlay->EndFrame(menuSettings, mousePosAtRightBtnDown, std::nullopt, spinnerVisible.load());
				mousePosAtRightBtnDown.reset();
				overlay->Render();
                glfwSwapBuffers(window);
            }
            continue;
        }
        if (!activeSimulationId || !renderContexts.contains(*activeSimulationId))
            activeSimulationId = renderContexts.begin()->first;
        activeRenderContext = &renderContexts.at(*activeSimulationId);
        ConsumeInputs();
        if (allowUserInputs) tiled = false;

        bool updatedPositions = false;
        bool anyNewTask = false;
        for (auto& [simulationId, currentRenderContext] : renderContexts) {
            auto& incomingRenderTasks = currentRenderContext.incomingRenderTasks;
            auto& currentRenderTask = currentRenderContext.currentRenderTask;
            const bool visible = tiled || simulationId == activeSimulationId;
            if (currentRenderContext.revolveCamera) {
                const auto now = std::chrono::high_resolution_clock::now();
                const float elapsedSeconds = std::chrono::duration<float>(now - currentRenderContext.lastRevolveTime).count();
                if (visible)
                    currentRenderContext.camera->Update(-elapsedSeconds * 2.f * PI / 5.f, 0.f, 0.f);
                currentRenderContext.lastRevolveTime = now;
            }

            bool newTask = false;
            if (!incomingRenderTasks.empty()) {
                auto incomingRenderTask = std::move(incomingRenderTasks.front());
                incomingRenderTasks.pop_front();
                if (const auto* update = std::get_if<std::unique_ptr<SimulationTaskUpdate>>(&incomingRenderTask)) {
                    if (*update && std::holds_alternative<std::unique_ptr<AtomRenderTask>>(currentRenderTask)) {
                        PrepareNewRenderTask(currentRenderContext, *std::get<std::unique_ptr<AtomRenderTask>>(currentRenderTask), **update);
                        updatedPositions = true;
                    }
                }
                else if (!std::holds_alternative<Rendering::NoTask>(incomingRenderTask)) {
                    if (dragSimulationId == simulationId) CancelInteraction();
                    currentRenderTask = std::move(incomingRenderTask);
                    newTask = true;
                }
            }
            if (newTask || currentRenderContext.shouldRecolorAtoms) {
                if (!std::holds_alternative<Rendering::NoTask>(currentRenderTask))
                    PrepareTask(currentRenderContext, currentRenderTask, !newTask);
                currentRenderContext.shouldRecolorAtoms = false;
                anyNewTask = true;
            }

			if (currentRenderContext.renderDataPipe
				&& std::holds_alternative<std::unique_ptr<AtomRenderTask>>(currentRenderTask)) {
				auto& atomTask = *std::get<std::unique_ptr<AtomRenderTask>>(currentRenderTask);
				currentRenderContext.renderPositionsHost.resize(currentRenderContext.renderDataPipe->PositionCount());
				int64_t renderStep = -1;
				if (currentRenderContext.renderDataPipe->TryCopyToHost(
					currentRenderContext.renderPositionsHost.data(), currentRenderContext.renderPositionsHost.size(), renderStep)) {
					auto status = atomTask.simStatus;
					status.step = renderStep;
					PrepareNewRenderTask(currentRenderContext, atomTask, Rendering::SimulationTaskUpdate{
						currentRenderContext.renderPositionsHost.data(), nullptr, status });
					updatedPositions = true;
				}
			}

        }

		std::vector<SimulationTab> tabs;
		tabs.reserve(renderContexts.size());
		for (auto& [simulationId, context] : renderContexts) {
			if (context.renderDataPipe) {
				bool completed = false;
				auto status = context.renderDataPipe->GetStatus(completed);
				context.completed = completed;
				if (std::holds_alternative<std::unique_ptr<AtomRenderTask>>(context.currentRenderTask))
					std::get<std::unique_ptr<AtomRenderTask>>(context.currentRenderTask)->simStatus = std::move(status);
			}
			tabs.push_back(SimulationTab{ simulationId,
				context.label.empty() ? std::format("Simulation {}", simulationId + 1) : context.label,
				simulationId == *activeSimulationId, context.completed });
		}



        const int msPerFrame = std::floor(1. / 60. * 1000.);
        bool shouldDraw = viewports.empty() || anyNewTask || updatedPositions || newInput || framebufferWasResized
			|| frameTime.elapsed().count() > msPerFrame;

        if (shouldDraw && framebufferSize.x > 0 && framebufferSize.y > 0) {
            RenderFrame(tabs);

			fps->NewFrame();
            frameTime = TimeIt{};
        }
    }
    activeRenderContext = nullptr;
    renderContexts.clear(); // Release GL resources while the render thread owns the GL context.
}

bool Display::RemoveStoppedRenderContexts() {
	bool removed = false;
	for (auto it = renderContexts.begin(); it != renderContexts.end();) {
		if (!it->second.renderDataPipe || it->second.renderDataPipe->GetState() != RenderDataPipe::State::Stopped) {
			++it;
			continue;
		}
		if (dragSimulationId == it->first)
			CancelInteraction();
		if (popupSimulationId == it->first)
			popupSimulationId.reset();
		if (activeSimulationId == it->first) {
			activeSimulationId.reset();
			activeRenderContext = nullptr;
		}
		viewports.erase(it->first);
		it = renderContexts.erase(it);
		removed = true;
	}
	return removed;
}

bool Display::ApplyPendingFramebufferResize() {
	if (!framebufferResizePending || framebufferSize.x <= 0 || framebufferSize.y <= 0)
		return false;

	CancelInteraction();
	viewports.clear();
	framebufferResizePending = false;
	if (activeRenderContext)
		activeRenderContext->camera->UpdateViewport(framebufferSize);
	glViewport(0, 0, framebufferSize.x, framebufferSize.y);
	if (renderTargetControl)
		renderTargetControl->Resize(framebufferSize);
	return true;
}

void Display::Submit(SimulationId simId, Rendering::Task task, bool blocking,
	std::shared_ptr<RenderDataPipe> renderDataPipe, std::string label) {
    {
        if (std::holds_alternative<std::unique_ptr<Rendering::SimulationTaskUpdate>>(task)) {
            if (std::get<std::unique_ptr<SimulationTaskUpdate>>(task) == nullptr) {
                int a = 0;
            }
        }
		std::lock_guard<std::mutex> lock(incomingRenderTaskMutex);

		incomingRenderTasksGlobal.push_back({simId, std::move(task), renderDataPipe, std::move(label)});
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
void Display::Free(SimulationId simId) {
	Submit(simId, Rendering::FreeTask{}, false);
}

void Display::UpdateSelection(SimulationId simId, const std::set<int>& selection,
	std::optional<Rendering::MoleculeInfo> selectedMolecule) {

	std::lock_guard<std::mutex> lock(inputMutex);
    newSelectionInputs.emplace_back(simId, selection, std::move(selectedMolecule));
}








bool Display::initGLFW() {
    // Initialize the library
    if (!glfwInit()) {
        throw std::runtime_error("\nGLFW failed to initialize");
    }

    GLFWmonitor* primaryMonitor = glfwGetPrimaryMonitor();
    if (!primaryMonitor) {
        glfwTerminate();
        return false;
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
	if (gizmoEnabled && activeRenderContext && activeRenderContext->activeGizmo
		&& (activeRenderContext->activeGizmo->pullForce || activeRenderContext->activeGizmo->rotateForce)) {
		auto& gizmo = *activeRenderContext->activeGizmo;
        return LiveEdit::MoveMolecule(Convert(gizmo.pullForce.value_or(glm::vec3{})), Convert(gizmo.rotateForce.value_or(glm::vec3{})));
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


	display.Submit(0, std::make_unique<Rendering::AtomRenderTask>(pclusters, pcMetas, params), true);
	display.Submit(0, std::make_unique<Rendering::SimulationTaskUpdate>(position.get(), nullptr, SimStatus{}), true);
}
