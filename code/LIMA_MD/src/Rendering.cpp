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

void Display::_RenderAtoms(Float3 boxSize, int totalParticles, bool fromCuda) {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	{
        const glm::mat4 VP = camera.ViewProjection();
        drawBoxOutlineShader->Draw(VP, Float3{ boxSize });
	}

    // TODO: Add coloringmethod flag, and let shaders discard a fragment if not showing solvents! (or just pass atomLetter colors as a buffer, where solvents can have alpha=0)
    const glm::mat4 view = camera.View();
    const glm::mat4 projection = camera.Projection();
    
    if (fromCuda)
        drawAtomsFromCudaShader->Draw(view, projection, totalParticles);
    else
        drawAtomsFromCpuShader->Draw(view, projection, totalParticles);	
}

void Display::PrepareNewRenderTask(const Rendering::SimulationTask& task)
{
    camera.Update(task.boxparams.BoxSizeFloat());

    auto start = std::chrono::high_resolution_clock::now();

    if (!drawBoxOutlineShader)
        drawBoxOutlineShader = std::make_unique<DrawBoxOutlineShader>();

    if (!drawAtomsFromCpuShader)
        drawAtomsFromCpuShader = std::make_unique<DrawAtomsShader<false>>(task.boxparams.totalParticles, nullptr, windowSize);


    //std::string windowText = window_title + "\n" + task.siminfo;
    //glfwSetWindowTitle(window, windowText.c_str());

    // Preprocess the renderAtoms
    {
        renderAtomsTemp.resize(task.boxparams.totalParticles);

        int index = 0;
        for (int pcid = 0; pcid < task.pcMeta.size(); pcid++) {
            for (int pid = 0; pid < 4; pid++) {
				const PersistentClusterMeta& pcMeta = task.pcMeta[pcid];

                if (pcMeta.particleIdsGlobal[pid] == -1)
					continue;

                auto atomType = RenderUtilities::RAS_getTypeFromAtomletter(pcMeta.atomLetter[pid]);
				const float chargeNormalized = (task.pclusters[pcid].pqd[pid].params.charge + elementaryChargeToKiloCoulombPerMole) / (elementaryChargeToKiloCoulombPerMole * 2.f); // I... think this might be bullshit/wrong?? :D
                renderAtomsTemp[index].position = task.positions[pcid * PersistentCluster::nParticles + pid].Tofloat4(RenderUtilities::getRadius(atomType));

                if (task.coloringMethod == ColoringMethod::Atomname)
                    renderAtomsTemp[index].color = RenderUtilities::getColor(atomType);
                else if (task.coloringMethod == ColoringMethod::Charge) {
                    renderAtomsTemp[index].color = RenderUtilities::GetColorInGradientBlueRed(chargeNormalized);
                }
                else if (task.coloringMethod == ColoringMethod::GradientFromCompoundId) {
                    renderAtomsTemp[index].color = RenderUtilities::GetColorInGradientHue(static_cast<float>(pcid) / task.pcMeta.size());
                }
                index++;
            }
        }
    }

    // Move the renderAtoms to device
    drawAtomsFromCpuShader->renderAtomsBuffer.SetData(renderAtomsTemp);
}



void Display::PrepareNewRenderTask(const Rendering::MoleculehullTask& task) {
    if (!drawBoxOutlineShader)
        drawBoxOutlineShader = std::make_unique<DrawBoxOutlineShader>();

    if (!drawTrianglesShader)
        drawTrianglesShader = std::make_unique<DrawTrianglesShader>();

    if (!drawAtomsFromCudaShader || drawAtomsFromCudaShader->numAtomsReservedInRenderatomsBuffer < task.molCollection.nParticles)
        drawAtomsFromCudaShader = std::make_unique<DrawAtomsShader<true>>(task.molCollection.nParticles, &renderAtomsBufferCudaResource, windowSize);

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
        drawAtomsFromCudaShader->Draw(V, P, molCollection.nParticles);

	if (renderFacets)
		drawTrianglesShader->Draw(VP, molCollection.facets, molCollection.nFacets, FacetDrawMode::EDGES, boxSize);

	if (renderFacetsNormals)
		drawNormalsShader->Draw(VP, molCollection.facets, molCollection.nFacets, boxSize);


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

	if (!drawAtomsFromCpuShader)
		drawAtomsFromCpuShader = std::make_unique<DrawAtomsShader<false>>(nAtoms, &renderAtomsBufferCudaResource, windowSize);



    camera.Update(task.grofile.box_size);

	// Preprocess the renderAtoms
	{
		renderAtomsTemp.resize(nAtoms);

        for (int i = 0; i < nAtoms; i++) {
			renderAtomsTemp[i].position = task.grofile.atoms[i].position.Tofloat4(RenderUtilities::getRadius(RenderUtilities::RAS_getTypeFromAtomletter(task.grofile.atoms[i].atomName[0])));

            if (task.highlightedAtoms.contains(i)) 
				renderAtomsTemp[i].color = float4(227.f / 255.f, 28.f / 255.f, 121.f / 255.f, 1.f) ; // Highlighted atoms are pink
            else if (task.coloringMethod == GradientFromAtomid)
                renderAtomsTemp[i].color = RenderUtilities::GetColorInGradientBlueRed(static_cast<float>(i) / nAtoms);
            else 
			    renderAtomsTemp[i].color = RenderUtilities::getColor(RenderUtilities::RAS_getTypeFromAtomletter(task.grofile.atoms[i].atomName[0]));
		}
	}

	// Move the renderAtoms to device
	{
        drawAtomsFromCpuShader->renderAtomsBuffer.SetData(renderAtomsTemp);
	}

}