#include "DisplayV2.h"
#include "Shaders.h"    

#include <GL/glew.h>
#include <GLFW/glfw3.h>


const float deg2rad = 2.f * PI / 360.f;
const float rad2deg = 1.f / deg2rad;

static glm::mat4 GetMVPMatrix(float camera_distance, float camera_pitch, float camera_yaw, int screenWidth, int screenHeight) {
    // Set up the perspective projection
    double aspectRatio = static_cast<double>(screenWidth) / static_cast<double>(screenHeight);
    double fovY = 45.0;
    double nearPlane = 0.1;
    double farPlane = 10.0;
    double fH = tan(fovY / 360.0 * 3.141592653589793) * nearPlane;
    double fW = fH * aspectRatio;

    glm::mat4 projection = glm::frustum(-fW, fW, -fH, fH, nearPlane, farPlane);

    // Set up the view matrix
    glm::mat4 view = glm::mat4(1.0f);
    view = glm::translate(view, glm::vec3(0.0f, 0.0f, camera_distance));
    view = glm::rotate(view, glm::radians(-90.0f), glm::vec3(1.0f, 0.0f, 0.0f));  // Fixed rotation around x-axis
    view = glm::rotate(view, glm::radians(camera_pitch), glm::vec3(1.0f, 0.0f, 0.0f));  // Rotation around x-axis for pitch
    view = glm::rotate(view, glm::radians(camera_yaw), glm::vec3(0.0f, 0.0f, 1.0f));  // Rotation around z-axis for yaw (y-axis in GLM's coordinate system)

    // Model matrix (identity in this case)
    glm::mat4 model = glm::mat4(1.0f);

    return projection * view * model;
}


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
                }                              
            }
        }
        };

    glfwSetKeyCallback(window, keyCallback);

    logger.finishSection("Display initialized");
}

Display::~Display() {
    glfwTerminate();
}

struct Circle {
    float2 center;
    float radius;
    float4 color; 
};

// Function to draw circles
void DrawBackgroundCircles() {
    Float3 baseColor{ 0x0e / 255.f, 0x1e / 255.f, 0x42 / 255.f };
    baseColor = baseColor * 0.89f;
    float2 center{ -1.f, 0.9f };


    glDisable(GL_DEPTH_TEST);  // Disable depth testing
    for (int i = 0; i < 5; i++) {
        Float3 color = baseColor * static_cast<float>(i) / 5.f;
        float4 colorf4 = color.Tofloat4(1.f);
        const float radius = 2.9f - 0.5f * static_cast<float>(i);
        glColor4f(colorf4.x, colorf4.y, colorf4.z, colorf4.w); // Set color

        glBegin(GL_TRIANGLE_FAN); // Use GL_TRIANGLE_FAN for circle drawing
        glVertex2f(center.x, center.y); // Center vertex

        int numSegments = 100; // Adjust for smoothness
        for (int i = 0; i <= numSegments; ++i) {
            float angle = 2.0f * PI * i / numSegments;
            float x = radius * cosf(angle) + center.x;
            float y = radius * sinf(angle) + center.y;
            glVertex2f(x, y);
        }

        glEnd();
    }
    glEnable(GL_DEPTH_TEST);  // Disable depth testing
}


void Display::render(const Float3* positions, const std::vector<Compound>& compounds, const BoxParams& boxparams, int64_t step, float temperature, ColoringMethod coloringMethod) {
    auto start = std::chrono::high_resolution_clock::now();

    if (!drawBoxOutlineShader)
        drawBoxOutlineShader = std::make_unique<DrawBoxOutlineShader>();

    if (!drawAtomsShader)
		drawAtomsShader = std::make_unique<DrawAtomsShader>(boxparams.total_particles, &renderAtomsBufferCudaResource);


    // Preprocess the renderAtoms
    {
        // Map buffer object for writing from CUDA
        RenderAtom* renderAtomsBuffer;
        cudaGraphicsMapResources(1, &renderAtomsBufferCudaResource, 0);
        size_t num_bytes = 0;
        cudaGraphicsResourceGetMappedPointer((void**)&renderAtomsBuffer, &num_bytes, renderAtomsBufferCudaResource);
        assert(num_bytes == boxparams.total_particles * sizeof(RenderAtom));

        rasterizer.render(positions, compounds, boxparams, step, camera_normal, coloringMethod, renderAtomsBuffer);

        // Release buffer object from CUDA
        cudaGraphicsUnmapResources(1, &renderAtomsBufferCudaResource, 0);
    }


    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


    //DrawBackgroundCircles();

    const glm::mat4 MVP = GetMVPMatrix(camera_distance, camera_pitch * rad2deg, camera_yaw * rad2deg, screenWidth, screenHeight);
    drawBoxOutlineShader->Draw(MVP);

    drawAtomsShader->Draw(MVP, boxparams.total_particles);

    // Swap front and back buffers
    glfwSwapBuffers(window);

    std::string window_text = std::format("{}        Step: {}    Temperature: {:.1f}[k]", window_title, step, temperature);
    glfwSetWindowTitle(window, window_text.c_str());

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);



    //printf("\tRender time: %4d ys  ", (int) duration.count());
}










void Display::Render(const MoleculeHullCollection& molCollection, Float3 boxSize) {
    if (!drawBoxOutlineShader)
        drawBoxOutlineShader = std::make_unique<DrawBoxOutlineShader>();

    if (!drawTrianglesShader)
        drawTrianglesShader = std::make_unique<DrawTrianglesShader>();

    if (!drawAtomsShader)
        drawAtomsShader = std::make_unique<DrawAtomsShader>(molCollection.nParticles, &renderAtomsBufferCudaResource);

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    const glm::mat4 MVP = GetMVPMatrix(camera_distance, camera_pitch * rad2deg, camera_yaw * rad2deg, screenWidth, screenHeight);
    drawBoxOutlineShader->Draw(MVP);

    {
        // Map buffer object for writing from CUDA
        RenderAtom* renderAtomsBuffer;
        cudaGraphicsMapResources(1, &renderAtomsBufferCudaResource, 0);
        size_t num_bytes = 0;
        cudaGraphicsResourceGetMappedPointer((void**)&renderAtomsBuffer, &num_bytes, renderAtomsBufferCudaResource);
        assert(num_bytes == molCollection.nParticles * sizeof(RenderAtom));

        cudaMemcpy(renderAtomsBuffer, molCollection.particles, sizeof(RenderAtom) * molCollection.nParticles, cudaMemcpyDeviceToDevice);

        // Release buffer object from CUDA
        cudaGraphicsUnmapResources(1, &renderAtomsBufferCudaResource, 0);
    }

    drawAtomsShader->Draw(MVP, molCollection.nParticles);

    drawTrianglesShader->Draw(MVP, molCollection.facets, molCollection.nFacets, FacetDrawMode::FACES, boxSize);


    // Swap front and back buffers
    glfwSwapBuffers(window);
}


void Display::Render(const std::vector<Facet>& facets, std::optional<Float3> facetsColor, bool drawFacetNormals, FacetDrawMode drawMode,
    const std::vector<Float3>& points, std::optional<Float3> pointsColor
    , Float3 boxSize, bool flushGLBuffer, bool swapGLBuffers) 
{
    if (!drawBoxOutlineShader)
        drawBoxOutlineShader = std::make_unique<DrawBoxOutlineShader>();

    if (!drawTrianglesShader)
        drawTrianglesShader = std::make_unique<DrawTrianglesShader>();

    if (!drawAtomsShader || drawAtomsShader->numAtoms < points.size())
        drawAtomsShader = std::make_unique<DrawAtomsShader>(points.size(), &renderAtomsBufferCudaResource);

    if (!drawNormalsShader)
        drawNormalsShader = std::make_unique<DrawNormalsShader>();

    if (flushGLBuffer)
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    const glm::mat4 MVP = GetMVPMatrix(camera_distance, camera_pitch * rad2deg, camera_yaw * rad2deg, screenWidth, screenHeight);
    drawBoxOutlineShader->Draw(MVP);

    if (points.size()) {
        // Map buffer object for writing from CUDA
        RenderAtom* renderAtomsBuffer;
        cudaGraphicsMapResources(1, &renderAtomsBufferCudaResource, 0);
        size_t num_bytes = 0;
        cudaGraphicsResourceGetMappedPointer((void**)&renderAtomsBuffer, &num_bytes, renderAtomsBufferCudaResource);
        assert(num_bytes >= points.size() * sizeof(RenderAtom));


        std::vector<RenderAtom> renderAtoms(points.size());
        for (int i = 0; i < points.size(); i++) {
            renderAtoms[i] = RenderAtom{ points[i], boxSize, 'O' };
            if (pointsColor.has_value())
                renderAtoms[i].color = pointsColor.value().Tofloat4(1.f);
		}

        cudaMemcpy(renderAtomsBuffer, renderAtoms.data(), sizeof(RenderAtom) * points.size(), cudaMemcpyHostToDevice);

        // Release buffer object from CUDA
        cudaGraphicsUnmapResources(1, &renderAtomsBufferCudaResource, 0);
    }

    drawAtomsShader->Draw(MVP, points.size());


    Facet* facetsDev;
    cudaMalloc(&facetsDev, sizeof(Facet) * facets.size());
    cudaMemcpy(facetsDev, facets.data(), sizeof(Facet) * facets.size(), cudaMemcpyHostToDevice);
    drawTrianglesShader->Draw(MVP, facetsDev, facets.size(), drawMode, boxSize);
    if (drawFacetNormals)
        drawNormalsShader->Draw(MVP, facetsDev, facets.size(), boxSize);
    cudaFree(facetsDev);


    // Swap front and back buffers
    if (swapGLBuffers)
        glfwSwapBuffers(window);
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

