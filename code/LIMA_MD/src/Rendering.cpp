#include "DisplayV2.h"
#include "Shaders.h"
#include "TimeIt.h"

//#include <GL/glew.h>
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





void Display::RenderLoop(const MoleculeHullCollection& molCollection, Float3 boxSize) {
    if (!drawBoxOutlineShader)
        drawBoxOutlineShader = std::make_unique<DrawBoxOutlineShader>();

    if (!drawTrianglesShader)
        drawTrianglesShader = std::make_unique<DrawTrianglesShader>();

    if (!drawAtomsShader)
        drawAtomsShader = std::make_unique<DrawAtomsShader>(molCollection.nParticles, &renderAtomsBufferCudaResource);




    while (true) {

        if (!checkWindowStatus())
            break;

        if (debugValue == 1) {
            debugValue = 0;
            break;
        }


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

        fps.NewFrame();
        std::string windowText = window_title + "    FPS: " + std::to_string(fps.GetFps());
        glfwSetWindowTitle(window, windowText.c_str());
    }
}


void CopyBufferIntoCudaIntoOpengl(cudaGraphicsResource** renderAtomsBufferCudaResource, Float3 boxSize, std::optional<Float3> pointsColor, const std::vector<Float3>& points) {
    // Map buffer object for writing from CUDA
    RenderAtom* renderAtomsBuffer;
    cudaGraphicsMapResources(1, renderAtomsBufferCudaResource, 0);
    size_t num_bytes = 0;
    cudaGraphicsResourceGetMappedPointer((void**)&renderAtomsBuffer, &num_bytes, *renderAtomsBufferCudaResource);
    assert(num_bytes >= points.size() * sizeof(RenderAtom));


    std::vector<RenderAtom> renderAtoms(points.size());
    for (int i = 0; i < points.size(); i++) {
        renderAtoms[i] = RenderAtom{ points[i], boxSize, 'O' };
        if (pointsColor.has_value())
            renderAtoms[i].color = pointsColor.value().Tofloat4(1.f);
    }

    cudaMemcpy(renderAtomsBuffer, renderAtoms.data(), sizeof(RenderAtom) * points.size(), cudaMemcpyHostToDevice);

    // Release buffer object from CUDA
    cudaGraphicsUnmapResources(1, renderAtomsBufferCudaResource, 0);
}


void Display::RenderLoop(std::vector<FacetTask>& facetTasks, std::vector<PointsTask>& pointsTasks, Float3 boxSize, std::optional<std::chrono::milliseconds> duration)
{
    if (!drawBoxOutlineShader)
        drawBoxOutlineShader = std::make_unique<DrawBoxOutlineShader>();

    if (!drawTrianglesShader)
        drawTrianglesShader = std::make_unique<DrawTrianglesShader>();

    int maxPoints = 0;
    for (const auto& task : pointsTasks)
		maxPoints = std::max(maxPoints, (int)task.points.size());

    if (!drawAtomsShader || drawAtomsShader->numAtoms < maxPoints)
        drawAtomsShader = std::make_unique<DrawAtomsShader>(maxPoints, &renderAtomsBufferCudaResource);

    if (!drawNormalsShader)
        drawNormalsShader = std::make_unique<DrawNormalsShader>();


    TimeIt timer{};

    while (true) {

        if (!checkWindowStatus())
            break;

        if (debugValue == 1) {
			debugValue = 0;
			break;
		}

        if (duration.has_value() && timer.elapsed() > duration)
			break;


        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        const glm::mat4 MVP = GetMVPMatrix(camera_distance, camera_pitch * rad2deg, camera_yaw * rad2deg, screenWidth, screenHeight);
        drawBoxOutlineShader->Draw(MVP);


        for (const auto& task : pointsTasks) {
            CopyBufferIntoCudaIntoOpengl(&renderAtomsBufferCudaResource, boxSize, task.pointsColor, task.points);
            drawAtomsShader->Draw(MVP, task.points.size());
        }

        for (const auto& task : facetTasks) {
            Facet* facetsDev;
            cudaMalloc(&facetsDev, sizeof(Facet) * task.facets.size());
            cudaMemcpy(facetsDev, task.facets.data(), sizeof(Facet) * task.facets.size(), cudaMemcpyHostToDevice);
            drawTrianglesShader->Draw(MVP, facetsDev, task.facets.size(), task.facetDrawMode, boxSize);
            if (task.drawFacetNormals)
                drawNormalsShader->Draw(MVP, facetsDev, task.facets.size(), boxSize);
            cudaFree(facetsDev);
        }


        // Swap front and back buffers
        glfwSwapBuffers(window);

        fps.NewFrame();
        std::string windowText = window_title + "    FPS: " + std::to_string(fps.GetFps());
        glfwSetWindowTitle(window, windowText.c_str());
    }
}