#include "DisplayV2.h"
#include "Shaders.h"
#include "TimeIt.h"

#include "RenderUtilities.cuh"
//#include <GL/glew.h>
#include <GLFW/glfw3.h>

const float deg2rad = 2.f * PI / 360.f;
const float rad2deg = 1.f / deg2rad;



static glm::mat4 GetMVPMatrix(float camera_distance, float camera_pitch, float camera_yaw, int screenWidth, int screenHeight, std::optional<Float3> boxSize = std::nullopt) {
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
    view = glm::rotate(view, glm::radians(camera_yaw), glm::vec3(0.0f, 0.0f, 1.0f));  // Rotation around z-axis for yaw

    // Model matrix (identity by default)
    glm::mat4 model = glm::mat4(1.0f);

    // Apply normalization transform if requested
    if (boxSize.has_value()) {
        //glm::vec3 boxCenter = boxSize.value().ToVec3() * 0.5f;
        glm::vec3 scale = 1.0f / boxSize.value().ToVec3();  // Normalizing scale

        // Translation to center the box at origin and scaling to normalize
        model = glm::translate(model, glm::vec3{ -0.5f }) * glm::scale(glm::mat4(1.0f), scale);
    }

    return projection * view * model;
}





void Display::_Render(const BoxParams& boxparams) {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	{
		const glm::mat4 MVP = GetMVPMatrix(camera_distance, camera_pitch * rad2deg, camera_yaw * rad2deg, screenWidth, screenHeight);
		drawBoxOutlineShader->Draw(MVP);
	}

	const glm::mat4 MVP = GetMVPMatrix(camera_distance, camera_pitch * rad2deg, camera_yaw * rad2deg, screenWidth, screenHeight, boxparams.boxSize);
	drawAtomsShader->Draw(MVP, boxparams.total_particles);

	glfwSwapBuffers(window);
}



void Display::PrepareNewRenderTask(const Float3* positions, const std::vector<Compound> compounds,
    const BoxParams boxparams, int64_t step, float temperature, ColoringMethod coloringMethod)
{
    if (boxparams.n_compounds < 0 || boxparams.n_compounds > 1000000)
        throw std::runtime_error("Invalid number of compounds");

    auto start = std::chrono::high_resolution_clock::now();

    if (!drawBoxOutlineShader)
        drawBoxOutlineShader = std::make_unique<DrawBoxOutlineShader>();

    if (!drawAtomsShader)
        drawAtomsShader = std::make_unique<DrawAtomsShader>(boxparams.total_particles, &renderAtomsBufferCudaResource);

    std::string window_text = std::format("{}        Step: {}    Temperature: {:.1f}[k]", window_title, step, temperature);
    glfwSetWindowTitle(window, window_text.c_str());

    // Preprocess the renderAtoms
    {
        renderAtomsTemp.resize(boxparams.total_particles);

        int index = 0;
        for (int cid = 0; cid < compounds.size(); cid++) {
            for (int pid = 0; pid < compounds[cid].n_particles; pid++) {
                auto atomType = RenderUtilities::RAS_getTypeFromAtomletter(compounds[cid].atomLetters[pid]);
                renderAtomsTemp[index].position = positions[cid * MAX_COMPOUND_PARTICLES + pid].Tofloat4(RenderUtilities::getRadius(atomType));

                if (coloringMethod == ColoringMethod::Atomname)
                    renderAtomsTemp[index].color = RenderUtilities::getColor(atomType);
                else if (coloringMethod == ColoringMethod::Charge) {
                    const float chargeNormalized = (static_cast<float>(compounds[cid].atom_charges[pid]) + elementaryChargeToKiloCoulombPerMole) / (elementaryChargeToKiloCoulombPerMole * 2.f);
                    renderAtomsTemp[index].color = float4{ chargeNormalized, 0.f, (1.f - chargeNormalized), 1.f };
                }

                /*      if (atomType == RenderUtilities::ATOM_TYPE::H && !drawHydrogens) {
                          renderAtomsTemp[index].position.w = 0.f;
                      }*/

                index++;
            }
        }
        for (int sid = 0; sid < boxparams.n_solvents; sid++) {
            renderAtomsTemp[index].position = positions[MAX_COMPOUND_PARTICLES * compounds.size() + sid].Tofloat4(RenderUtilities::getRadius(RenderUtilities::ATOM_TYPE::SOL));

            renderAtomsTemp[index].color = RenderUtilities::getColor(RenderUtilities::ATOM_TYPE::SOL);
            index++;
        }
    }

    // Move the renderAtoms to device
    {
        // Map buffer object for writing from CUDA
        RenderAtom* renderAtomsBuffer;
        cudaGraphicsMapResources(1, &renderAtomsBufferCudaResource, 0);
        size_t num_bytes = 0;
        cudaGraphicsResourceGetMappedPointer((void**)&renderAtomsBuffer, &num_bytes, renderAtomsBufferCudaResource);

        if (num_bytes != boxparams.total_particles * sizeof(RenderAtom)) {
            throw std::runtime_error("RenderAtom buffer size mismatch");
        }

        assert(num_bytes == boxparams.total_particles * sizeof(RenderAtom));

        //rasterizer->render(positions, boxparams, step, coloringMethod, renderAtomsBuffer);

        cudaMemcpy(renderAtomsBuffer, renderAtomsTemp.data(), sizeof(RenderAtom) * renderAtomsTemp.size(), cudaMemcpyHostToDevice);

        // Release buffer object from CUDA
        cudaGraphicsUnmapResources(1, &renderAtomsBufferCudaResource, 0);
    }
}



void Display::PrepareNewRenderTask(const MoleculeHullCollection& molCollection) {
    if (!drawBoxOutlineShader)
        drawBoxOutlineShader = std::make_unique<DrawBoxOutlineShader>();

    if (!drawTrianglesShader)
        drawTrianglesShader = std::make_unique<DrawTrianglesShader>();

    if (!drawAtomsShader || drawAtomsShader->numAtomsReservedInRenderatomsBuffer < molCollection.nParticles)
        drawAtomsShader = std::make_unique<DrawAtomsShader>(molCollection.nParticles, &renderAtomsBufferCudaResource);

    if (!drawNormalsShader)
        drawNormalsShader = std::make_unique<DrawNormalsShader>();

    if (renderAtoms) {
        // Map buffer object for writing from CUDA
        RenderAtom* renderAtomsBuffer;
        cudaGraphicsMapResources(1, &renderAtomsBufferCudaResource, 0);
        size_t num_bytes = 0;

        cudaGraphicsResourceGetMappedPointer((void**)&renderAtomsBuffer, &num_bytes, renderAtomsBufferCudaResource);
        assert(num_bytes >= molCollection.nParticles * sizeof(RenderAtom));

        cudaMemcpy(renderAtomsBuffer, molCollection.particles, sizeof(RenderAtom) * molCollection.nParticles, cudaMemcpyDeviceToDevice);

        // Release buffer object from CUDA
        cudaGraphicsUnmapResources(1, &renderAtomsBufferCudaResource, 0);
    }
}



void Display::_Render(const MoleculeHullCollection& molCollection, Float3 boxSize) {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	const glm::mat4 MVP = GetMVPMatrix(camera_distance, camera_pitch * rad2deg, camera_yaw * rad2deg, screenWidth, screenHeight);
	drawBoxOutlineShader->Draw(MVP);

    if (renderAtoms) 
        drawAtomsShader->Draw(MVP, molCollection.nParticles);

	if (renderFacets)
		drawTrianglesShader->Draw(MVP, molCollection.facets, molCollection.nFacets, FacetDrawMode::EDGES, boxSize);

	if (renderFacetsNormals)
		drawNormalsShader->Draw(MVP, molCollection.facets, molCollection.nFacets, boxSize);

	// Swap front and back buffers
	glfwSwapBuffers(window);

	fps.NewFrame();
	std::string windowText = window_title + "    FPS: " + std::to_string(fps.GetFps());
	glfwSetWindowTitle(window, windowText.c_str());
}


//void CopyBufferIntoCudaIntoOpengl(cudaGraphicsResource** renderAtomsBufferCudaResource, Float3 boxSize, std::optional<Float3> pointsColor, const std::vector<Float3>& points) {
//    // Map buffer object for writing from CUDA
//    RenderAtom* renderAtomsBuffer;
//    cudaGraphicsMapResources(1, renderAtomsBufferCudaResource, 0);
//    size_t num_bytes = 0;
//    cudaGraphicsResourceGetMappedPointer((void**)&renderAtomsBuffer, &num_bytes, *renderAtomsBufferCudaResource);
//    assert(num_bytes >= points.size() * sizeof(RenderAtom));
//
//
//    std::vector<RenderAtom> renderAtoms(points.size());
//    for (int i = 0; i < points.size(); i++) {
//        renderAtoms[i] = RenderAtom{ points[i], boxSize, 'O' };
//        if (pointsColor.has_value())
//            renderAtoms[i].color = pointsColor.value().Tofloat4(1.f);
//    }
//
//    cudaMemcpy(renderAtomsBuffer, renderAtoms.data(), sizeof(RenderAtom) * points.size(), cudaMemcpyHostToDevice);
//
//    // Release buffer object from CUDA
//    cudaGraphicsUnmapResources(1, renderAtomsBufferCudaResource, 0);
//}

//
//void Display::RenderLoop(std::vector<FacetTask>& facetTasks, std::vector<PointsTask>& pointsTasks, Float3 boxSize, std::optional<std::chrono::milliseconds> duration)
//{
//    if (!drawBoxOutlineShader)
//        drawBoxOutlineShader = std::make_unique<DrawBoxOutlineShader>();
//
//    if (!drawTrianglesShader)
//        drawTrianglesShader = std::make_unique<DrawTrianglesShader>();
//
//    int maxPoints = 0;
//    for (const auto& task : pointsTasks)
//		maxPoints = std::max(maxPoints, (int)task.points.size());
//
//    if (!drawAtomsShader || drawAtomsShader->numAtomsReservedInRenderatomsBuffer < maxPoints)
//        drawAtomsShader = std::make_unique<DrawAtomsShader>(maxPoints, &renderAtomsBufferCudaResource);
//
//    if (!drawNormalsShader)
//        drawNormalsShader = std::make_unique<DrawNormalsShader>();
//
//
//    TimeIt timer{};
//
//    while (true) {
//
//        //if (!checkWindowStatus())
//        //    break;
//
//        if (debugValue == 1) {
//			debugValue = 0;
//			break;
//		}
//
//        if (duration.has_value() && timer.elapsed() > duration)
//			break;
//
//
//        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
//
//        const glm::mat4 MVP = GetMVPMatrix(camera_distance, camera_pitch * rad2deg, camera_yaw * rad2deg, screenWidth, screenHeight);
//        drawBoxOutlineShader->Draw(MVP);
//
//
//        for (const auto& task : pointsTasks) {
//            CopyBufferIntoCudaIntoOpengl(&renderAtomsBufferCudaResource, boxSize, task.pointsColor, task.points);
//            drawAtomsShader->Draw(MVP, task.points.size());
//        }
//
//        for (const auto& task : facetTasks) {
//            Facet* facetsDev;
//            cudaMalloc(&facetsDev, sizeof(Facet) * task.facets.size());
//            cudaMemcpy(facetsDev, task.facets.data(), sizeof(Facet) * task.facets.size(), cudaMemcpyHostToDevice);
//            drawTrianglesShader->Draw(MVP, facetsDev, task.facets.size(), task.facetDrawMode, boxSize);
//            if (task.drawFacetNormals)
//                drawNormalsShader->Draw(MVP, facetsDev, task.facets.size(), boxSize);
//            cudaFree(facetsDev);
//        }
//
//
//        // Swap front and back buffers
//        glfwSwapBuffers(window);
//
//        fps.NewFrame();
//        std::string windowText = window_title + "    FPS: " + std::to_string(fps.GetFps());
//        glfwSetWindowTitle(window, windowText.c_str());
//    }
//}