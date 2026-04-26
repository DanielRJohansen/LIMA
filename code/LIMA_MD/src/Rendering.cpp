#include "Display.h"
#include "Shaders.h"
#include "TimeIt.h"
#include "MDFiles.h"

#include "RenderUtilities.cuh"
//#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include "SSBO.h"

const float deg2rad = 2.f * PI / 360.f;
const float rad2deg = 1.f / deg2rad;


glm::mat4 Camera::View() const {
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

glm::mat4 Camera::Projection() const {
	//double aspectRatio = static_cast<double>(screenWidth) / static_cast<double>(screenHeight);
	double aspectRatio = 1.f;
	double fovY = 45.0;
	double nearPlane = 0.1;
	double farPlane = 1000.0;
	double fH = tan(glm::radians(fovY / 2.0)) * nearPlane;
	double fW = fH * aspectRatio;

	return glm::frustum(-fW, fW, -fH, fH, nearPlane, farPlane);
}

glm::mat4 Camera::ViewProjection() const {
	return Projection() * View();
}

glm::vec3 AnyPerpendicular(const glm::vec3& dir)
{
	const glm::vec3 helper = std::abs(dir.z) < 0.999f
		? glm::vec3(0.f, 0.f, 1.f)
		: glm::vec3(0.f, 1.f, 0.f);
	return glm::normalize(glm::cross(helper, dir));
}

glm::mat4 RotationFromZAxisTo(const glm::vec3& direction)
{
	const glm::vec3 localAxis(0.f, 0.f, 1.f);
	const float c = glm::clamp(glm::dot(localAxis, direction), -1.f, 1.f);

	if (c > 0.99999f)
		return glm::mat4(1.f);

	if (c < -0.99999f) {
		const glm::vec3 rotAxis = AnyPerpendicular(localAxis);
		return glm::rotate(glm::mat4(1.f), PI, rotAxis);
	}

	const glm::vec3 rotAxis = glm::normalize(glm::cross(localAxis, direction));
	const float angle = std::acos(c);
	return glm::rotate(glm::mat4(1.f), angle, rotAxis);
}

Arrow::Arrow(glm::vec3 direction, glm::vec4 color, int id) : direction(glm::normalize(direction)), color(color), uniqueId(id) {
	constexpr int radialSegments = 64;
	constexpr float totalLength = 1.0f;
	constexpr float shaftLength = 0.78f * totalLength;
	constexpr float shaftRadius = 0.035f * totalLength;
	constexpr float headLength = totalLength - shaftLength;
	constexpr float headRadius = 0.09f * totalLength;

	const glm::vec3 axis(0.f, 0.f, 1.f);
	const glm::vec3 u(1.f, 0.f, 0.f);
	const glm::vec3 v(0.f, 1.f, 0.f);

	const glm::vec3 shaftStart(0.f, 0.f, 0.f);
	const glm::vec3 shaftEnd(0.f, 0.f, shaftLength);
	const glm::vec3 coneBase(0.f, 0.f, shaftLength);
	const glm::vec3 apex(0.f, 0.f, totalLength);

	vertices.reserve(radialSegments * 12);

	auto AppendTriangle = [&](const glm::vec3& a, const glm::vec3& b, const glm::vec3& c) {
		const glm::vec3 normal = glm::normalize(glm::cross(b - a, c - a));
		vertices.push_back({ a, normal });
		vertices.push_back({ b, normal });
		vertices.push_back({ c, normal });
		};

	for (int i = 0; i < radialSegments; ++i) {
		const float a0 = 2.f * 3.1415f * static_cast<float>(i) / static_cast<float>(radialSegments);
		const float a1 = 2.f * 3.1415f * static_cast<float>(i + 1) / static_cast<float>(radialSegments);

		const glm::vec3 r0s = std::cos(a0) * u * shaftRadius + std::sin(a0) * v * shaftRadius;
		const glm::vec3 r1s = std::cos(a1) * u * shaftRadius + std::sin(a1) * v * shaftRadius;

		const glm::vec3 p0 = shaftStart + r0s;
		const glm::vec3 p1 = shaftStart + r1s;
		const glm::vec3 q0 = shaftEnd + r0s;
		const glm::vec3 q1 = shaftEnd + r1s;

		AppendTriangle(p0, p1, q1);
		AppendTriangle(p0, q1, q0);

		const glm::vec3 r0c = std::cos(a0) * u * headRadius + std::sin(a0) * v * headRadius;
		const glm::vec3 r1c = std::cos(a1) * u * headRadius + std::sin(a1) * v * headRadius;

		const glm::vec3 c0 = coneBase + r0c;
		const glm::vec3 c1 = coneBase + r1c;

		AppendTriangle(c0, c1, apex);
		AppendTriangle(coneBase, c1, c0);
	}

	AppendTriangle(shaftStart, shaftStart + glm::vec3(shaftRadius, 0.f, 0.f), shaftStart + glm::vec3(0.f, shaftRadius, 0.f));
}


void Arrow::Draw(DrawTrianglesShader* shader, const glm::mat4& VP, const glm::vec3& pos, float scale) const {
	//const float length = 2.f;
	const glm::vec3 localAxis(0.f, 0.f, 1.f);

	glm::mat4 rotation(1.f);
	const float c = glm::clamp(glm::dot(localAxis, direction), -1.f, 1.f);

	if (c < 0.99999f) {
		if (c > -0.99999f) {
			const glm::vec3 rotAxis = glm::normalize(glm::cross(localAxis, direction));
			const float angle = std::acos(c);
			rotation = glm::rotate(glm::mat4(1.f), angle, rotAxis);
		}
		else {
			const glm::vec3 rotAxis = AnyPerpendicular(localAxis);
			rotation = glm::rotate(glm::mat4(1.f), 3.1415f, rotAxis);
		}
	}

	glm::mat4 model(1.f);
	model = glm::translate(model, pos);
	model *= rotation;
	model = glm::scale(model, glm::vec3(scale));

	glm::mat4 MVP = VP * model;

	shader->Draw(vertices, MVP, model, color, uniqueId);
};

Ring::Ring(glm::vec3 normal, glm::vec4 color, int id)
	: normal(glm::normalize(normal)), color(color), uniqueId(id)
{
	constexpr int majorSegments = 96;
	constexpr int minorSegments = 12;
	constexpr float majorRadius = 1.45f;
	constexpr float tubeRadius = 0.020f;

	vertices.reserve(majorSegments * minorSegments * 6);

	auto TorusPoint = [&](float u, float v) {
		const float cu = std::cos(u);
		const float su = std::sin(u);
		const float cv = std::cos(v);
		const float sv = std::sin(v);

		const float r = majorRadius + tubeRadius * cv;
		return glm::vec3(r * cu, r * su, tubeRadius * sv);
		};

	auto TorusNormal = [&](float u, float v) {
		const float cu = std::cos(u);
		const float su = std::sin(u);
		const float cv = std::cos(v);
		const float sv = std::sin(v);

		return glm::normalize(glm::vec3(cv * cu, cv * su, sv));
		};

	auto AppendTri = [&](const glm::vec3& a, const glm::vec3& na,
		const glm::vec3& b, const glm::vec3& nb,
		const glm::vec3& c, const glm::vec3& nc)
		{
			vertices.push_back({ a, na });
			vertices.push_back({ b, nb });
			vertices.push_back({ c, nc });
		};

	for (int i = 0; i < majorSegments; ++i) {
		const float u0 = 2.f * PI * static_cast<float>(i) / static_cast<float>(majorSegments);
		const float u1 = 2.f * PI * static_cast<float>(i + 1) / static_cast<float>(majorSegments);

		for (int j = 0; j < minorSegments; ++j) {
			const float v0 = 2.f * PI * static_cast<float>(j) / static_cast<float>(minorSegments);
			const float v1 = 2.f * PI * static_cast<float>(j + 1) / static_cast<float>(minorSegments);

			const glm::vec3 p00 = TorusPoint(u0, v0);
			const glm::vec3 p10 = TorusPoint(u1, v0);
			const glm::vec3 p11 = TorusPoint(u1, v1);
			const glm::vec3 p01 = TorusPoint(u0, v1);

			const glm::vec3 n00 = TorusNormal(u0, v0);
			const glm::vec3 n10 = TorusNormal(u1, v0);
			const glm::vec3 n11 = TorusNormal(u1, v1);
			const glm::vec3 n01 = TorusNormal(u0, v1);

			AppendTri(p00, n00, p10, n10, p11, n11);
			AppendTri(p00, n00, p11, n11, p01, n01);
		}
	}
}

void Ring::Draw(DrawTrianglesShader* shader, const glm::mat4& VP, const glm::vec3& pos, float scale) const
{
	glm::mat4 model(1.f);
	model = glm::translate(model, pos);
	model *= RotationFromZAxisTo(normal);
	model = glm::scale(model, glm::vec3(scale));

	const glm::mat4 MVP = VP * model;
	shader->Draw(vertices, MVP, model, color, uniqueId);
}

void TransformGizmo::Draw(DrawTrianglesShader* shader, const glm::mat4& VP) const {
	const int axis = activeAxis.value_or(-1);

	const float translateScaleX = activeMode == GizmoMode::Translate && axis == 0 ? 2.2f : 2.f;
	const float translateScaleY = activeMode == GizmoMode::Translate && axis == 1 ? 2.2f : 2.f;
	const float translateScaleZ = activeMode == GizmoMode::Translate && axis == 2 ? 2.2f : 2.f;

	const float rotateScaleX = activeMode == GizmoMode::Rotate && axis == 0 ? 2.2f : 2.f;
	const float rotateScaleY = activeMode == GizmoMode::Rotate && axis == 1 ? 2.2f : 2.f;
	const float rotateScaleZ = activeMode == GizmoMode::Rotate && axis == 2 ? 2.2f : 2.f;

	arrowX.Draw(shader, VP, position, translateScaleX);
	arrowY.Draw(shader, VP, position, translateScaleY);
	arrowZ.Draw(shader, VP, position, translateScaleZ);

	ringX.Draw(shader, VP, position, rotateScaleX);
	ringY.Draw(shader, VP, position, rotateScaleY);
	ringZ.Draw(shader, VP, position, rotateScaleZ);
}



















void Display::_RenderAtoms() {
	
	const glm::mat4 VP = camera.ViewProjection();

	// TODO: Add coloringmethod flag, and let shaders discard a fragment if not showing solvents! (or just pass atomLetter colors as a buffer, where solvents can have alpha=0)
	const glm::mat4 view = camera.View();
	const glm::mat4 projection = camera.Projection();

	//drawAtomsPrettyShader->Draw(*renderAtomsBuffer, renderAtomsHost.size(), view, projection);	
	drawAtomsFromCpuShader->Draw(*renderAtomsBuffer, renderAtomsHost.size(), view, projection);
}

int Display::GetObjectIdAtPixel(glm::ivec2 pixel)
{
	auto scopedDrawBinding = renderTargetControl->BindForDraw();
	renderTargetControl->ClearForPicking();

	_RenderAtoms();

	// Must be done last!
	if (activeGizmo) {
		glClear(GL_DEPTH_BUFFER_BIT);   // forget scene depth
		activeGizmo->Draw(drawTrianglesShader.get(), camera.ViewProjection());
	}
	int elementId = renderTargetControl->ReadIdAtPixel(pixel);
	//printf("ElementId %d\n", elementId);
	return elementId;
}

void Display::PrepareNewRenderTask(const Rendering::SimulationTask& task, bool ignorePosition)
{
	camera.Update(task.boxparams.BoxSizeFloat());

	if (!drawBoxOutlineShader)
		drawBoxOutlineShader = std::make_unique<DrawBoxOutlineShader>();
	if (!drawAtomsFromCpuShader)
		drawAtomsFromCpuShader = std::make_unique<DrawAtomsShader>();
	if (!drawTrianglesShader)
		drawTrianglesShader = std::make_unique<DrawTrianglesShader>();
	if (!renderTargetControl)
		renderTargetControl = std::make_unique<RenderTargetControl>();
	if (!drawAtomsPrettyShader)
		drawAtomsPrettyShader = std::make_unique<DrawAtomsPrettyShader>();
	renderTargetControl->Resize(windowSize);

	// Preprocess the renderAtoms
	{
		renderAtomsHost.resize(task.boxparams.totalParticles, RenderAtom{});
		for (int pcid = 0; pcid < task.pcMeta.size(); pcid++) {
			for (int pid = 0; pid < 4; pid++) {
				const PersistentClusterMeta& pcMeta = task.pcMeta[pcid];
				const int pidGlobal = pcMeta.particleIdsGlobal[pid];

				if (pidGlobal == -1)
					continue;

				auto atomType = RenderUtilities::RAS_getTypeFromAtomletter(pcMeta.atomLetter[pid], pcMeta.isSolvent);
				const float chargeNormalized = (task.pclusters[pcid].pqd[pid].params.charge + elementaryChargeToKiloCoulombPerMole) / (elementaryChargeToKiloCoulombPerMole * 2.f); // I... think this might be bullshit/wrong?? :D

				if (!ignorePosition)
					renderAtomsHost[pidGlobal].position = task.pclusters[pcid].pqd[pid].position.Tofloat4(RenderUtilities::getRadius(atomType));
				renderAtomsHost[pidGlobal].flags.y = pcMeta.particleIdsGlobal[pid];

				if (rendersettings.coloringMethod == ColoringMethod::Atomname)
					renderAtomsHost[pidGlobal].color = RenderUtilities::getColor(atomType);
				else if (rendersettings.coloringMethod == ColoringMethod::Charge) {
					renderAtomsHost[pidGlobal].color = RenderUtilities::GetColorInGradientBlueRed(chargeNormalized);
				}
				else if (rendersettings.coloringMethod == ColoringMethod::PersistentClusterId) {
					int nElementsPerRevolution = 12;
					float fraction = (static_cast<float>(pcid % nElementsPerRevolution) / static_cast<float>(nElementsPerRevolution));
					renderAtomsHost[pidGlobal].color = RenderUtilities::GetColorInGradientHue(fraction);
				}
				else if (rendersettings.coloringMethod == ColoringMethod::GradientFromAtomid) {					
					renderAtomsHost[pidGlobal].color = RenderUtilities::GetColorInGradientHue(static_cast<float>(pidGlobal) / static_cast<float>(task.boxparams.totalParticles));
				}

				if (!rendersettings.showSolvents && pcMeta.isSolvent)
					renderAtomsHost[pidGlobal].color.w = 0.f;
			}
		}
	}

	if (activeGizmo && activeGizmo->idOfAtomAttachedTo != -1 && activeGizmo->idOfAtomAttachedTo < renderAtomsHost.size()) {
		int attachedAtomId = activeGizmo->idOfAtomAttachedTo;
		if (attachedAtomId < renderAtomsHost.size()) {
			activeGizmo->position = glm::vec3(renderAtomsHost[attachedAtomId].position.x, renderAtomsHost[attachedAtomId].position.y, renderAtomsHost[attachedAtomId].position.z);
		}
	}

	// Move the renderAtoms to device
	renderAtomsBuffer->SetData(renderAtomsHost);
}

void Display::PrepareNewRenderTask(Rendering::SimulationTask& currentTask, const Rendering::SimulationTaskUpdate& update)
{
	currentTask.simStatus = update.simStatus;

	// Update the renderAtoms
	{
		for (int pcid = 0; pcid < currentTask.pcMeta.size(); pcid++) {
			for (int pid = 0; pid < 4; pid++) {
				const PersistentClusterMeta& pcMeta = currentTask.pcMeta[pcid];
				const int pidGlobal = pcMeta.particleIdsGlobal[pid];
				if (pidGlobal == -1)
					continue;
				renderAtomsHost[pidGlobal].position = update.positions[pcid * PersistentCluster::maxParticles + pid].Tofloat4(renderAtomsHost[pidGlobal].position.w);
			}
		}
	}
	if (activeGizmo && activeGizmo->idOfAtomAttachedTo != -1 && activeGizmo->idOfAtomAttachedTo < renderAtomsHost.size()) {
		int attachedAtomId = activeGizmo->idOfAtomAttachedTo;
		if (attachedAtomId < renderAtomsHost.size()) {
			activeGizmo->position = glm::vec3(renderAtomsHost[attachedAtomId].position.x, renderAtomsHost[attachedAtomId].position.y, renderAtomsHost[attachedAtomId].position.z);
		}
	}
	// Move the renderAtoms to device
	renderAtomsBuffer->SetData(renderAtomsHost);
}



void Display::PrepareNewRenderTask(const Rendering::MoleculehullTask& task) {
	if (!drawBoxOutlineShader)
		drawBoxOutlineShader = std::make_unique<DrawBoxOutlineShader>();

	if (!drawFacetsShader)
		drawFacetsShader = std::make_unique<DrawFacetsShader>();

	/*if (!drawAtomsFromCudaShader || drawAtomsFromCudaShader->numAtomsReservedInRenderatomsBuffer < task.molCollection.nParticles)
		drawAtomsFromCudaShader = std::make_unique<DrawAtomsShader<true>>(task.molCollection.nParticles, &renderAtomsBufferCudaResource, windowSize);*/

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
	//const glm::mat4 MVP = GetMVPMatrix(camera_distance, camera_pitch * rad2deg, camera_yaw * rad2deg, screenWidth, screenHeight, boxSize.x);
	const glm::mat4 V = camera.View();
	const glm::mat4 P = camera.Projection();
	const glm::mat4 VP = camera.ViewProjection();

	/*if (renderAtoms)
		drawAtomsFromCudaShader->Draw(*renderAtomsBuffer, V, P, molCollection.nParticles);*/

	if (renderFacets)
		drawFacetsShader->Draw(VP, molCollection.facets, molCollection.nFacets, FacetDrawMode::EDGES, boxSize);

	if (renderFacetsNormals)
		drawNormalsShader->Draw(VP, molCollection.facets, molCollection.nFacets, boxSize);


	fps.NewFrame();
	std::string windowText = window_title + "    FPS: " + std::to_string(fps.GetFps());
	glfwSetWindowTitle(window, windowText.c_str());
}

void Display::_Render(const Rendering::Task& currentRenderTask) {

	// Check shaders is Init
	if (!drawBackgroundGradientShader)
		drawBackgroundGradientShader = std::make_unique<DrawBackgroundGradientShader>();
	if (!drawBoxOutlineShader)
		drawBoxOutlineShader = std::make_unique<DrawBoxOutlineShader>();

	SimStatus simStatus{};
	Float3 boxSize{};

	// First extract necessary information from the render task
	if (!std::holds_alternative<Rendering::NoTask>(currentRenderTask)) {
		std::visit([&](auto& taskPtr) {
			using T = std::decay_t<decltype(taskPtr)>;
			if constexpr (std::is_same_v<T, std::unique_ptr<Rendering::SimulationTask>>) {
				const int nParticles = taskPtr->boxparams.totalParticles;
				simStatus = taskPtr->simStatus;
				boxSize = taskPtr->boxparams.BoxSizeFloat();
				//_RenderAtoms(taskPtr->boxparams.BoxSizeFloat(), nParticles, false);
			}
			else if constexpr (std::is_same_v<T, std::unique_ptr<Rendering::MoleculehullTask>>) {
				//_Render(taskPtr->molCollection, taskPtr->boxSize);
				boxSize = taskPtr->boxSize;
			}
			else if constexpr (std::is_same_v<T, std::unique_ptr<Rendering::GrofileTask>>) {
				//_RenderAtoms(taskPtr->grofile.box_size, taskPtr->nAtoms, false);
				boxSize = taskPtr->grofile.box_size;
			}
			}, currentRenderTask);
	}

	const glm::mat4 VP = camera.ViewProjection();


	// START OF RENDERING
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	drawBackgroundGradientShader->Draw(ColorScheme::backgroundBot, ColorScheme::backgroundTop);
	drawBoxOutlineShader->Draw(VP, boxSize);
	

	// Then render
	if (!std::holds_alternative<Rendering::NoTask>(currentRenderTask)) {
		std::visit([&](auto& taskPtr) {
			using T = std::decay_t<decltype(taskPtr)>;
			if constexpr (std::is_same_v<T, std::unique_ptr<Rendering::SimulationTask>>) {
				const int nParticles = taskPtr->boxparams.totalParticles;
				_RenderAtoms();
			}
			else if constexpr (std::is_same_v<T, std::unique_ptr<Rendering::MoleculehullTask>>) {
				_Render(taskPtr->molCollection, taskPtr->boxSize);
			}
			else if constexpr (std::is_same_v<T, std::unique_ptr<Rendering::GrofileTask>>) {
				_RenderAtoms();
			}
			}, currentRenderTask);
	}

	
	if (activeGizmo) {
		// DO NOT RENDER ANYTHING IN 3D AFTER THIS POINT
		glClear(GL_DEPTH_BUFFER_BIT);   // forget scene depth
		activeGizmo->Draw(drawTrianglesShader.get(), VP);
		glEnable(GL_DEPTH_TEST);
	}

	overlay->enableConsole = allowUserInputs;
	overlay->Draw(rendersettings, simStatus, fps.GetFps(), mousePosAtRightBtnDown);
	mousePosAtRightBtnDown = std::nullopt;
	overlay->Render();

	glfwSwapBuffers(window);
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
		drawAtomsFromCpuShader = std::make_unique<DrawAtomsShader>();



	camera.Update(task.grofile.box_size);

	// Preprocess the renderAtoms
	{
		renderAtomsHost.resize(nAtoms);

		for (int i = 0; i < nAtoms; i++) {
			renderAtomsHost[i].position = task.grofile.atoms[i].position.Tofloat4(RenderUtilities::getRadius(RenderUtilities::RAS_getTypeFromAtomletter(task.grofile.atoms[i].atomName[0])));

			if (task.highlightedAtoms.contains(i))
				renderAtomsHost[i].color = float4(227.f / 255.f, 28.f / 255.f, 121.f / 255.f, 1.f); // Highlighted atoms are pink
			else if (rendersettings.coloringMethod == ColoringMethod::GradientFromAtomid)
				renderAtomsHost[i].color = RenderUtilities::GetColorInGradientBlueRed(static_cast<float>(i) / nAtoms);
			else
				renderAtomsHost[i].color = RenderUtilities::getColor(RenderUtilities::RAS_getTypeFromAtomletter(task.grofile.atoms[i].atomName[0]));
		}
	}

	// Move the renderAtoms to device
	{
		renderAtomsBuffer->SetData(renderAtomsHost);
	}

}

void Display::_UpdateSelection(const std::set<int>& selection) {
	// This is purposefully done in 2 passes, as the selection is likely MUCH smaller that the renderatoms, and this no point in doing lookings.
	for (auto& atom : renderAtomsHost) {
		atom.HighLight(false);
	}
	for (int id : selection) {
		if (id < renderAtomsHost.size())
			renderAtomsHost[id].HighLight(true);
	}
	renderAtomsBuffer->SetData(renderAtomsHost);
}
