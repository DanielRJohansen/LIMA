#include "DisplayV2.h"
#include <glfw3.h>


#ifndef __linux__



void Display::updateCamera(float pitch, float yaw, float delta_distance) {

    // Reset previous rotations
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    camera_distance += delta_distance;
    glTranslatef(0.0f, .0f, camera_distance);  // Move the scene away from the camera    
    glRotatef(-90.0f, 1.0f, 0.0f, 0.0f);  // Rotate around the x-axis
    

    // Apply rotation for pitch and yaw
    // Rotate around the x-axis for pitch
    glRotatef(pitch*rad2deg, 1.0f, 0.0f, 0.0f);
    // Rotate around the y-axis for yaw
    glRotatef(yaw*rad2deg, 0.0f, 0.0f, 1.0f);
    
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

    glfwSetWindowUserPointer(window, this);

    auto keyCallback = [](GLFWwindow* window, int key, int scancode, int action, int mods) {
        if (action == GLFW_PRESS) {
            // Retrieve the Display instance from the window user pointer
            Display* display = static_cast<Display*>(glfwGetWindowUserPointer(window));
            if (display) {
                const float delta = 3.1415 / 4.f;
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
                }                              
            }
        }
        };

    glfwSetKeyCallback(window, keyCallback);

    logger.finishSection("Display initialized");
}

Display::~Display() {
    TerminateGLEW();
    glfwTerminate();
    TerminateGLEW();
}


void Display::render(const Float3* positions, const std::vector<Compound>& compounds, const BoxParams& boxparams, int64_t step, float temperature, ColoringMethod coloringMethod) {
    auto start = std::chrono::high_resolution_clock::now();

    if (!pipelineInitialized) {
        initializePipeline(boxparams.total_particles);
        pipelineInitialized = true;
	}



    // Preprocess the renderAtoms
    
    // Map buffer object for writing from CUDA
    RenderAtom* renderAtomsBuffer;
    cudaGraphicsMapResources(1, &renderAtomsBufferCudaResource, 0);
    size_t num_bytes = 0;
	cudaGraphicsResourceGetMappedPointer((void**)&renderAtomsBuffer, &num_bytes, renderAtomsBufferCudaResource);
    assert(num_bytes == boxparams.total_particles * sizeof(RenderAtom));

    rasterizer.render(positions, compounds, boxparams, step, camera_normal, coloringMethod, renderAtomsBuffer);

    // Release buffer object from CUDA
    cudaGraphicsUnmapResources(1, &renderAtomsBufferCudaResource, 0);
    

    

    glClear(GL_COLOR_BUFFER_BIT);


    DrawBoxOutline();
    DrawAtoms(boxparams.total_particles);

    // Swap front and back buffers
    glfwSwapBuffers(window);


    std::string window_text = std::format("{}        Step: {}    Temperature: {:.1f}[k]", window_title, step, temperature);
    glfwSetWindowTitle(window, window_text.c_str());

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    //printf("\tRender time: %4d ys  ", (int) duration.count());
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
        throw std::exception("\nGLFW failed to initialize");
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

#endif
