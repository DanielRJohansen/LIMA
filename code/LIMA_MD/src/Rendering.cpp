#include "Display.h"
#include "Shaders.h"
#include "TimeIt.h"
#include "MDFiles.h"

#include "RenderUtilities.cuh"
#include "NewCartoonRenderer.h"
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
	if (rendersettings.coloringMethod == ColoringMethod::NewCartoon
		&& newCartoonRenderer && newCartoonRenderer->HasGeometry()) {
		newCartoonRenderer->Draw(*drawTrianglesShader, VP);
		return;
	}

	// TODO: Add coloringmethod flag, and let shaders discard a fragment if not showing solvents! (or just pass atomLetter colors as a buffer, where solvents can have alpha=0)
	const glm::mat4 view = camera.View();
	const glm::mat4 projection = camera.Projection();

	//drawAtomsPrettyShader->Draw(*renderAtomsBuffer, renderAtomsHost.size(), view, projection);	
	drawAtomsFromCpuShader->Draw(*renderAtomsBuffer, renderAtomsHost.size(), view, projection);
}

int Display::GetObjectIdAtPixel(glm::ivec2 pixel)
{
	if (!renderTargetControl || windowSize.x <= 0 || windowSize.y <= 0
		|| framebufferSize.x <= 0 || framebufferSize.y <= 0)
		return -1;

	// GLFW cursor positions are logical window coordinates; the picking
	// attachment uses framebuffer pixels.
	pixel.x = static_cast<int>(static_cast<double>(pixel.x) * framebufferSize.x / windowSize.x);
	pixel.y = static_cast<int>(static_cast<double>(pixel.y) * framebufferSize.y / windowSize.y);
	pixel.x = std::clamp(pixel.x, 0, framebufferSize.x - 1);
	pixel.y = std::clamp(pixel.y, 0, framebufferSize.y - 1);

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

Rendering::AtomRenderTask::AtomRenderTask(const GroFile& grofile, bool shouldShowSolvents)
	: positions(grofile.atoms.size())
	, atoms(grofile.atoms.size())
	, packedPositionIndices(grofile.atoms.size())
	, boxSize(grofile.box_size)
	, showSolvents(shouldShowSolvents)
{
	for (std::size_t atomId = 0; atomId < grofile.atoms.size(); ++atomId) {
		const GroRecord& atom = grofile.atoms[atomId];
		positions[atomId] = atom.position;
		packedPositionIndices[atomId] = static_cast<int>(atomId);
		atoms[atomId].atomLetter = atom.atomName[0];
		atoms[atomId].isSolvent = atom.residueName == "SOL" || atom.residueName == "TIP3";
	}
}

Rendering::AtomRenderTask::AtomRenderTask(
	const std::vector<PersistentCluster>& pclusters,
	const std::vector<PersistentClusterMeta>& pcMeta,
	const BoxParams& boxparams,
	SimStatus initialSimStatus,
	BackboneChains initialBackboneChains)
	: positions(boxparams.totalParticles)
	, atoms(boxparams.totalParticles)
	, packedPositionIndices(boxparams.totalParticles, -1)
	, boxSize(boxparams.BoxSizeFloat())
	, simStatus(initialSimStatus)
	, backboneChains(std::move(initialBackboneChains))
{
	for (std::size_t pcid = 0; pcid < pcMeta.size(); ++pcid) {
		for (int pid = 0; pid < PersistentCluster::maxParticles; ++pid) {
			const int globalParticleId = pcMeta[pcid].particleIdsGlobal[pid];
			if (globalParticleId < 0)
				continue;

			const std::size_t atomId = static_cast<std::size_t>(globalParticleId);
			if (atomId >= atoms.size() || pcid >= pclusters.size())
				throw std::runtime_error("Invalid simulation particle mapping in render task");

			positions[atomId] = pclusters[pcid].pqd[pid].position;
			packedPositionIndices[atomId] = static_cast<int>(pcid * PersistentCluster::maxParticles + pid);
			atoms[atomId] = {
				pcMeta[pcid].atomLetter[pid],
				pclusters[pcid].pqd[pid].params.charge,
				static_cast<int>(pcid),
				pcMeta[pcid].isSolvent
			};
		}
	}
}

void Display::PrepareNewRenderTask(Rendering::AtomRenderTask& task, bool ignorePosition)
{
	if (!ignorePosition)
		rendersettings.showSolvents = task.showSolvents;

	if (rendersettings.coloringMethod == ColoringMethod::NewCartoon) {
		if (!newCartoonRenderer)
			newCartoonRenderer = std::make_unique<NewCartoon::Renderer>();
		newCartoonRenderer->Prepare(
			task.backboneChains, task.positions, task.boxSize);
	}
	else if (newCartoonRenderer) {
		newCartoonRenderer->Clear();
	}

	camera.Update(task.boxSize);

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
	renderTargetControl->Resize(framebufferSize);

	// Preprocess the renderAtoms
	{
		renderAtomsHost.resize(task.atoms.size(), RenderAtom{});
		for (std::size_t atomId = 0; atomId < task.atoms.size(); ++atomId) {
			const Rendering::AtomRenderData& atom = task.atoms[atomId];
			const auto atomType = RenderUtilities::RAS_getTypeFromAtomletter(atom.atomLetter, atom.isSolvent);
			const float chargeNormalized = (atom.charge + elementaryChargeToKiloCoulombPerMole) / (elementaryChargeToKiloCoulombPerMole * 2.f);

			if (!ignorePosition)
				renderAtomsHost[atomId].position = task.positions[atomId].Tofloat4(RenderUtilities::getRadius(atomType));
			renderAtomsHost[atomId].flags.y = static_cast<int>(atomId);

			if (task.highlightedAtoms.contains(static_cast<int>(atomId)))
				renderAtomsHost[atomId].color = float4(227.f / 255.f, 28.f / 255.f, 121.f / 255.f, 1.f);
			else if (rendersettings.coloringMethod == ColoringMethod::Atomname
				|| rendersettings.coloringMethod == ColoringMethod::NewCartoon)
				renderAtomsHost[atomId].color = RenderUtilities::getColor(atomType);
			else if (rendersettings.coloringMethod == ColoringMethod::Charge)
				renderAtomsHost[atomId].color = RenderUtilities::GetColorInGradientBlueRed(chargeNormalized);
			else if (rendersettings.coloringMethod == ColoringMethod::PersistentClusterId) {
				constexpr int nElementsPerRevolution = 12;
				const int groupId = std::max(atom.groupId, 0);
				const float fraction = static_cast<float>(groupId % nElementsPerRevolution) / nElementsPerRevolution;
				renderAtomsHost[atomId].color = RenderUtilities::GetColorInGradientHue(fraction);
			}
			else if (rendersettings.coloringMethod == ColoringMethod::GradientFromAtomid)
				renderAtomsHost[atomId].color = RenderUtilities::GetColorInGradientHue(static_cast<float>(atomId) / task.atoms.size());
			else if (rendersettings.coloringMethod == ColoringMethod::ForceMagnitude)
				renderAtomsHost[atomId].color = RenderUtilities::GetLogColorGradient(0, 1e3f, 1e6f);

			if (!rendersettings.showSolvents && atom.isSolvent)
				renderAtomsHost[atomId].color.w = 0.f;
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

void Display::PrepareNewRenderTask(Rendering::AtomRenderTask& currentTask, const Rendering::SimulationTaskUpdate& update)
{
	currentTask.simStatus = update.simStatus;

	// Update the renderAtoms
	{
		for (std::size_t atomId = 0; atomId < currentTask.atoms.size(); ++atomId) {
			const int packedPositionIndex = currentTask.packedPositionIndices[atomId];
			if (packedPositionIndex < 0)
				continue;
			currentTask.positions[atomId] = update.positions[packedPositionIndex];
			renderAtomsHost[atomId].position = currentTask.positions[atomId].Tofloat4(renderAtomsHost[atomId].position.w);
			if (update.forceMagnitudes && rendersettings.coloringMethod == ColoringMethod::ForceMagnitude)
				renderAtomsHost[atomId].color = RenderUtilities::GetLogColorGradient(update.forceMagnitudes[atomId], 1e5f, 1e11f);
		}
	}
	if (rendersettings.coloringMethod == ColoringMethod::NewCartoon && newCartoonRenderer)
		newCartoonRenderer->Update(currentTask.positions);
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
	if (newCartoonRenderer)
		newCartoonRenderer->Clear();
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
	glViewport(0, 0, framebufferSize.x, framebufferSize.y);

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
			if constexpr (std::is_same_v<T, std::unique_ptr<Rendering::AtomRenderTask>>) {
				simStatus = taskPtr->simStatus;
				boxSize = taskPtr->boxSize;
			}
			else if constexpr (std::is_same_v<T, std::unique_ptr<Rendering::MoleculehullTask>>) {
				//_Render(taskPtr->molCollection, taskPtr->boxSize);
				boxSize = taskPtr->boxSize;
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
			if constexpr (std::is_same_v<T, std::unique_ptr<Rendering::AtomRenderTask>>) {
				_RenderAtoms();
			}
			else if constexpr (std::is_same_v<T, std::unique_ptr<Rendering::MoleculehullTask>>) {
				_Render(taskPtr->molCollection, taskPtr->boxSize);
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
	overlay->Draw(rendersettings, simStatus, fps.GetFps(), mousePosAtRightBtnDown,
		spinnerVisible.load());
	mousePosAtRightBtnDown = std::nullopt;
	overlay->Render();

	glfwSwapBuffers(window);
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
