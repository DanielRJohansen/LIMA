#include "Display.h"
#include "Shaders.h"
#include "TimeIt.h"
#include "MDFiles.h"

#include "RenderUtilities.cuh"
//#include <GL/glew.h>
#include <GLFW/glfw3.h>

const float deg2rad = 2.f * PI / 360.f;
const float rad2deg = 1.f / deg2rad;


glm::mat4 Camera::View() {
    glm::mat4 view = glm::mat4(1.0f);

    // Translate the camera back by the camera distance
    view = glm::translate(view, glm::vec3(0.0f, 0.0f, dist));

    // Apply the fixed rotation to make Z up
    view = glm::rotate(view, (-PI / 2.f), glm::vec3(1.0f, 0.0f, 0.0f));

    // Apply pitch and yaw rotations
    view = glm::rotate(view, pitch, glm::vec3(1.0f, 0.0f, 0.0f));  // Rotation around x-axis for pitch
    view = glm::rotate(view, yaw, glm::vec3(0.0f, 0.0f, 1.0f));    // Rotation around z-axis for yaw

    // Translate the world to the opposite direction of the camera position to look at the center
    view = glm::translate(view, ToVec3(-center));

    return view;
}

glm::mat4 Camera::Projection() {
    //double aspectRatio = static_cast<double>(screenWidth) / static_cast<double>(screenHeight);
    double aspectRatio = 1.f;
    double fovY = 45.0;
    double nearPlane = 0.1;
    double farPlane = 1000.0;
    double fH = tan(glm::radians(fovY / 2.0)) * nearPlane;
    double fW = fH * aspectRatio;

    return glm::frustum(-fW, fW, -fH, fH, nearPlane, farPlane);
}

glm::mat4 Camera::ViewProjection() {
    return Projection() * View();
}

void Display::_RenderAtomsFromCudaresource(Float3 boxSize, int totalParticles) {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	{
        const glm::mat4 VP = camera.ViewProjection();
        drawBoxOutlineShader->Draw(VP, Float3{ boxSize });
	}

    const glm::mat4 view = camera.View();
    const glm::mat4 projection = camera.Projection();
    

    drawAtomsShader->Draw(view, projection, totalParticles);

	glfwSwapBuffers(window);
}



void Display::PrepareNewRenderTask(const Rendering::SimulationTask& task)
{
    camera.Update(Float3{ task.boxparams.boxSize });

    if (task.boxparams.n_compounds < 0 || task.boxparams.n_compounds > 1000000)
        throw std::runtime_error("Invalid number of compounds");

    auto start = std::chrono::high_resolution_clock::now();

    if (!drawBoxOutlineShader)
        drawBoxOutlineShader = std::make_unique<DrawBoxOutlineShader>();

    if (!drawAtomsShader)
        drawAtomsShader = std::make_unique<DrawAtomsShader>(task.boxparams.total_particles, &renderAtomsBufferCudaResource);


    std::string windowText = window_title + "\n" + task.siminfo;
    glfwSetWindowTitle(window, windowText.c_str());

    // Preprocess the renderAtoms
    {
        renderAtomsTemp.resize(task.boxparams.total_particles);

        int index = 0;
        for (int cid = 0; cid < task.compounds.size(); cid++) {
            for (int pid = 0; pid < task.compounds[cid].n_particles; pid++) {
                auto atomType = RenderUtilities::RAS_getTypeFromAtomletter(task.compounds[cid].atomLetters[pid]);
                renderAtomsTemp[index].position = task.positions[cid * MAX_COMPOUND_PARTICLES + pid].Tofloat4(RenderUtilities::getRadius(atomType));

                if (task.coloringMethod == ColoringMethod::Atomname)
                    renderAtomsTemp[index].color = RenderUtilities::getColor(atomType);
                else if (task.coloringMethod == ColoringMethod::Charge) {
                    const float chargeNormalized = (static_cast<float>(task.compounds[cid].atom_charges[pid]) + elementaryChargeToKiloCoulombPerMole) / (elementaryChargeToKiloCoulombPerMole * 2.f);
                    renderAtomsTemp[index].color = RenderUtilities::GetColorInGradientBlueRed(chargeNormalized);
                }
                else if (task.coloringMethod == ColoringMethod::GradientFromCompoundId) {
                    renderAtomsTemp[index].color = RenderUtilities::GetColorInGradientHue(static_cast<float>(cid) / task.boxparams.n_compounds);
                }
                index++;
            }
        }
        for (int sid = 0; sid < task.boxparams.n_solvents; sid++) {
            renderAtomsTemp[index].position = task.positions[MAX_COMPOUND_PARTICLES * task.compounds.size() + sid].Tofloat4(RenderUtilities::getRadius(RenderUtilities::ATOM_TYPE::SOL));

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

        if (num_bytes != task.boxparams.total_particles * sizeof(RenderAtom)) {
            throw std::runtime_error("RenderAtom buffer size mismatch");
        }

        assert(num_bytes == task.boxparams.total_particles * sizeof(RenderAtom));

        cudaMemcpy(renderAtomsBuffer, renderAtomsTemp.data(), sizeof(RenderAtom) * renderAtomsTemp.size(), cudaMemcpyHostToDevice);

        // Release buffer object from CUDA
        cudaGraphicsUnmapResources(1, &renderAtomsBufferCudaResource, 0);
    }
}



void Display::PrepareNewRenderTask(const Rendering::MoleculehullTask& task) {
    if (!drawBoxOutlineShader)
        drawBoxOutlineShader = std::make_unique<DrawBoxOutlineShader>();

    if (!drawTrianglesShader)
        drawTrianglesShader = std::make_unique<DrawTrianglesShader>();

    if (!drawAtomsShader || drawAtomsShader->numAtomsReservedInRenderatomsBuffer < task.molCollection.nParticles)
        drawAtomsShader = std::make_unique<DrawAtomsShader>(task.molCollection.nParticles, &renderAtomsBufferCudaResource);

    if (!drawNormalsShader)
        drawNormalsShader = std::make_unique<DrawNormalsShader>();

    camera.Update(task.boxSize);

    if (renderAtoms) {
        // Map buffer object for writing from CUDA
        RenderAtom* renderAtomsBuffer;
        cudaGraphicsMapResources(1, &renderAtomsBufferCudaResource, 0);
        size_t num_bytes = 0;

        cudaGraphicsResourceGetMappedPointer((void**)&renderAtomsBuffer, &num_bytes, renderAtomsBufferCudaResource);
        assert(num_bytes >= task.molCollection.nParticles * sizeof(RenderAtom));

        cudaMemcpy(renderAtomsBuffer, task.molCollection.particles, sizeof(RenderAtom) * task.molCollection.nParticles, cudaMemcpyDeviceToDevice);

        // Release buffer object from CUDA
        cudaGraphicsUnmapResources(1, &renderAtomsBufferCudaResource, 0);
    }
}



void Display::_Render(const MoleculeHullCollection& molCollection, Float3 boxSize) {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	//const glm::mat4 MVP = GetMVPMatrix(camera_distance, camera_pitch * rad2deg, camera_yaw * rad2deg, screenWidth, screenHeight, boxSize.x);
    const glm::mat4 V = camera.View();
    const glm::mat4 P = camera.Projection();
    const glm::mat4 VP = camera.ViewProjection();
	drawBoxOutlineShader->Draw(VP, boxSize);

    if (renderAtoms) 
        drawAtomsShader->Draw(V, P, molCollection.nParticles);

	if (renderFacets)
		drawTrianglesShader->Draw(VP, molCollection.facets, molCollection.nFacets, FacetDrawMode::EDGES, boxSize);

	if (renderFacetsNormals)
		drawNormalsShader->Draw(VP, molCollection.facets, molCollection.nFacets, boxSize);

	// Swap front and back buffers
	glfwSwapBuffers(window);

	fps.NewFrame();
	std::string windowText = window_title + "    FPS: " + std::to_string(fps.GetFps());
	glfwSetWindowTitle(window, windowText.c_str());
}


void Display::PrepareNewRenderTask(Rendering::GrofileTask& task) {
    int nAtoms = task.grofile.atoms.size();
    if (!task.drawSolvent) {
        for (int i = 0; i < task.grofile.atoms.size(); i++) {
            auto resname = task.grofile.atoms[i].residueName;
            if (resname == "SOL" || resname == "TIP3") {
                nAtoms = i;
                break;
            }
        }
    }
    task.nAtoms = nAtoms;

	if (!drawBoxOutlineShader)
		drawBoxOutlineShader = std::make_unique<DrawBoxOutlineShader>();

	if (!drawAtomsShader)
		drawAtomsShader = std::make_unique<DrawAtomsShader>(nAtoms, &renderAtomsBufferCudaResource);



    camera.Update(task.grofile.box_size);

	// Preprocess the renderAtoms
	{
		renderAtomsTemp.resize(nAtoms);

        for (int i = 0; i < nAtoms; i++) {
			renderAtomsTemp[i].position = task.grofile.atoms[i].position.Tofloat4(RenderUtilities::getRadius(RenderUtilities::RAS_getTypeFromAtomletter(task.grofile.atoms[i].atomName[0])));

            if (task.coloringMethod == GradientFromAtomid)
                renderAtomsTemp[i].color = RenderUtilities::GetColorInGradientBlueRed(static_cast<float>(i) / nAtoms);
            else 
			    renderAtomsTemp[i].color = RenderUtilities::getColor(RenderUtilities::RAS_getTypeFromAtomletter(task.grofile.atoms[i].atomName[0]));
		}
	}

	// Move the renderAtoms to device
	{
		// Map buffer object for writing from CUDA
		RenderAtom* renderAtomsBuffer;
		cudaGraphicsMapResources(1, &renderAtomsBufferCudaResource, 0);
		size_t num_bytes = 0;
		cudaGraphicsResourceGetMappedPointer((void**)&renderAtomsBuffer, &num_bytes, renderAtomsBufferCudaResource);

		if (num_bytes != nAtoms * sizeof(RenderAtom)) {
			throw std::runtime_error("RenderAtom buffer size mismatch");
		}

		cudaMemcpy(renderAtomsBuffer, renderAtomsTemp.data(), sizeof(RenderAtom) * renderAtomsTemp.size(), cudaMemcpyHostToDevice);

		// Release buffer object from CUDA
		cudaGraphicsUnmapResources(1, &renderAtomsBufferCudaResource, 0);
	}

}