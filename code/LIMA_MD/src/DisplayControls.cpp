#include "Display.h"
#include "DisplayInternal.h"
#include "imgui.h"

#include <GLFW/glfw3.h>
#include <algorithm>
#include <ranges>




glm::vec3 GetAxisDirection(int axis)
{
    switch (axis) {
    case 0: return glm::vec3(1.f, 0.f, 0.f);
    case 1: return glm::vec3(0.f, 1.f, 0.f);
    case 2: return glm::vec3(0.f, 0.f, 1.f);
    default: return glm::vec3(0.f);
    }
}
glm::vec3 GetCameraWorldPosition(const glm::mat4& view)
{
    const glm::mat4 invView = glm::inverse(view);
    return glm::vec3(invView[3]);
}

glm::vec3 ScreenToWorldRayDirection(
    glm::vec2 mousePos,
    const glm::mat4& view,
    const glm::mat4& proj,
    glm::vec2 windowSize)
{
    const float x = 2.f * mousePos.x / windowSize.x - 1.f;
    const float y = 1.f - 2.f * mousePos.y / windowSize.y;

    const glm::vec4 rayClip{ x, y, -1.f, 1.f };

    glm::vec4 rayEye = glm::inverse(proj) * rayClip;
    rayEye = glm::vec4(rayEye.x, rayEye.y, -1.f, 0.f);

    return glm::normalize(glm::vec3(glm::inverse(view) * rayEye));
}

bool ClosestPointBetweenLines(
    const glm::vec3& p1, const glm::vec3& d1,
    const glm::vec3& p2, const glm::vec3& d2,
    float& t1Out)
{
    const glm::vec3 r = p1 - p2;
    const float a = glm::dot(d1, d1);
    const float e = glm::dot(d2, d2);
    const float b = glm::dot(d1, d2);
    const float c = glm::dot(d1, r);
    const float f = glm::dot(d2, r);

    const float denom = a * e - b * b;
    if (std::abs(denom) < 1e-8f)
        return false;

    t1Out = (b * f - c * e) / denom;
    return true;
}

bool IntersectRayPlane(
    const glm::vec3& rayOrigin,
    const glm::vec3& rayDir,
    const glm::vec3& planePoint,
    const glm::vec3& planeNormal,
    glm::vec3& hitPoint)
{
    const float denom = glm::dot(rayDir, planeNormal);
    if (std::abs(denom) < 1e-6f)
        return false;

    const float t = glm::dot(planePoint - rayOrigin, planeNormal) / denom;
    if (t < 0.f)
        return false;

    hitPoint = rayOrigin + rayDir * t;
    return true;
}

glm::vec3 ProjectOntoPlane(const glm::vec3& v, const glm::vec3& normal)
{
    return v - glm::dot(v, normal) * normal;
}

void TransformGizmo::SetActiveAxis(int selectedObjectId)
{
    if (selectedObjectId == arrowX.uniqueId) {
        activeAxis = 0;
        activeMode = GizmoMode::Translate;
    }
    else if (selectedObjectId == arrowY.uniqueId) {
        activeAxis = 1;
        activeMode = GizmoMode::Translate;
    }
    else if (selectedObjectId == arrowZ.uniqueId) {
        activeAxis = 2;
        activeMode = GizmoMode::Translate;
    }
    else if (selectedObjectId == ringX.uniqueId) {
        activeAxis = 0;
        activeMode = GizmoMode::Rotate;
    }
    else if (selectedObjectId == ringY.uniqueId) {
        activeAxis = 1;
        activeMode = GizmoMode::Rotate;
    }
    else if (selectedObjectId == ringZ.uniqueId) {
        activeAxis = 2;
        activeMode = GizmoMode::Rotate;
    }
    else {
        activeAxis = std::nullopt;
    }
}

void TransformGizmo::BeginDragging(glm::vec2 mousePos, const Camera& camera, glm::vec2 windowSize)
{
    if (!activeAxis)
        return;

    dragStartPosition = position;
    dragStartMousePos = mousePos;
    pullForce = std::nullopt;
    rotateForce = std::nullopt;
    dragStartRotateVector = glm::vec3(0.f);

    if (activeMode == GizmoMode::Translate) {
        pullForce = glm::vec3(0.f);
        return;
    }

    const glm::vec3 axisDir = GetAxisDirection(*activeAxis);
    const glm::mat4 view = camera.View();
    const glm::mat4 proj = camera.Projection();

    const glm::vec3 rayOrigin = GetCameraWorldPosition(view);
    const glm::vec3 rayDir = ScreenToWorldRayDirection(mousePos, view, proj, windowSize);

    glm::vec3 hitPoint{};
    if (!IntersectRayPlane(rayOrigin, rayDir, position, axisDir, hitPoint)) {
        rotateForce = glm::vec3(0.f);
        return;
    }

    glm::vec3 v = hitPoint - position;
    v = ProjectOntoPlane(v, axisDir);
    if (glm::length(v) < 1e-4f) {
        rotateForce = glm::vec3(0.f);
        return;
    }

    dragStartRotateVector = glm::normalize(v);
    rotateForce = glm::vec3(0.f);
}

void TransformGizmo::UpdateDraggingForce(glm::vec2 mousePos, const Camera& camera, glm::vec2 windowSize)
{
    if (!activeAxis.has_value()) {
        pullForce = std::nullopt;
        rotateForce = std::nullopt;
        return;
    }

    const glm::vec3 axisDir = GetAxisDirection(*activeAxis);

    if (activeMode == GizmoMode::Translate) {
        rotateForce = std::nullopt;

        const glm::vec3 axisOrigin = dragStartPosition;
        const glm::mat4 view = camera.View();
        const glm::mat4 proj = camera.Projection();

        const glm::vec3 rayOrigin = GetCameraWorldPosition(view);
        const glm::vec3 rayDirStart = ScreenToWorldRayDirection(dragStartMousePos, view, proj, windowSize);
        const glm::vec3 rayDirCurrent = ScreenToWorldRayDirection(mousePos, view, proj, windowSize);

        float axisTAtStart = 0.f;
        float axisTNow = 0.f;

        const bool okStart = ClosestPointBetweenLines(axisOrigin, axisDir, rayOrigin, rayDirStart, axisTAtStart);
        const bool okNow = ClosestPointBetweenLines(axisOrigin, axisDir, rayOrigin, rayDirCurrent, axisTNow);

        if (!okStart || !okNow) {
            pullForce = glm::vec3(0.f);
            return;
        }

        const float deltaAxis = axisTNow - axisTAtStart;
        const glm::vec3 targetPosition = dragStartPosition + axisDir * deltaAxis;

        const float axisError = glm::dot(targetPosition - position, axisDir);

        const float stiffness = 1.0f;
        float scalarForce = axisError * stiffness;

        const float maxForce = 1.0f;
        scalarForce = std::clamp(scalarForce, -maxForce, maxForce);

        pullForce = axisDir * scalarForce;
        return;
    }

    pullForce = std::nullopt;

    const glm::mat4 view = camera.View();
    const glm::mat4 proj = camera.Projection();

    const glm::vec3 rayOrigin = GetCameraWorldPosition(view);
    const glm::vec3 rayDir = ScreenToWorldRayDirection(mousePos, view, proj, windowSize);

    glm::vec3 hitPoint{};
    if (!IntersectRayPlane(rayOrigin, rayDir, dragStartPosition, axisDir, hitPoint)) {
        rotateForce = glm::vec3(0.f);
        return;
    }

    glm::vec3 currentVector = hitPoint - dragStartPosition;
    currentVector = ProjectOntoPlane(currentVector, axisDir);

    if (glm::length(currentVector) < 1e-4f || glm::length(dragStartRotateVector) < 1e-4f) {
        rotateForce = glm::vec3(0.f);
        return;
    }

    currentVector = glm::normalize(currentVector);

    const float sinAngle = glm::dot(axisDir, glm::cross(dragStartRotateVector, currentVector));
    const float cosAngle = glm::clamp(glm::dot(dragStartRotateVector, currentVector), -1.f, 1.f);
    const float signedAngle = std::atan2(sinAngle, cosAngle);

    const float stiffness = 2.0f;
    float scalarForce = signedAngle * stiffness;

    const float maxForce = 1.0f;
    scalarForce = std::clamp(scalarForce, -maxForce, maxForce);

    rotateForce = axisDir * scalarForce;
}







// ----------------------------------------- GLFW callbacks ----------------------------------------- //
void Display::OnMouseMove(double xpos, double ypos) {
	if (ImGui::GetCurrentContext() && ImGui::GetIO().WantCaptureMouse) {
		mousePos = { xpos, ypos };
		return;
	}
	if (!activeRenderContext)
		return;
	auto& renderContext = *activeRenderContext;
	if (gizmoEnabled && renderContext.activeGizmo && renderContext.activeGizmo->activeAxis.has_value()) {
		renderContext.activeGizmo->UpdateDraggingForce(glm::vec2(xpos, ypos), *renderContext.camera, windowSize);
    }
    else if (isDragging) {
		renderContext.revolveCamera = false;
        const float sensitivity = 0.001f;
        const float xOffset = static_cast<float>(xpos - mousePos.x) * sensitivity;
        const float yOffset = static_cast<float>(mousePos.y - ypos) * sensitivity;

		renderContext.camera->Update(xOffset, -yOffset, 0);
    }
	mousePos = { xpos, ypos };
}

void Display::HandleGizmo(int objectId) {
	if (!allowUserInputs || !gizmoEnabled || !activeRenderContext)
        return;
	auto& renderContext = *activeRenderContext;

    if (objectId == -1) {
		renderContext.activeGizmo.reset();
        stopMovingLiveeditCmd.store(true);
        return;
    }

    if (!renderContext.activeGizmo) {
		renderContext.activeGizmo = std::make_unique<TransformGizmo>();
    }
    if (activeRenderContext && objectId >= 0 && objectId < activeRenderContext->renderAtomsHost.size()) {
		renderContext.activeGizmo->position = glm::vec3{ activeRenderContext->renderAtomsHost[objectId].position.x, activeRenderContext->renderAtomsHost[objectId].position.y, activeRenderContext->renderAtomsHost[objectId].position.z };
		renderContext.activeGizmo->idOfAtomAttachedTo = objectId;
    }
}

void Display::SelectMolecule(int atomId) {
	if (!activeRenderContext)
		return;
	auto& renderContext = *activeRenderContext;
	if (!std::holds_alternative<std::unique_ptr<Rendering::AtomRenderTask>>(renderContext.currentRenderTask))
		return;
	const auto& molecules = std::get<std::unique_ptr<Rendering::AtomRenderTask>>(renderContext.currentRenderTask)->molecules;
	const auto molecule = std::ranges::find_if(molecules, [atomId](const auto& atoms) {
		return std::ranges::find(atoms.atomIds, atomId) != atoms.atomIds.end();
	});
	if (molecule == molecules.end()) {
		renderContext.selectedMolecule.reset();
		_UpdateSelection(renderContext, {});
		return;
	}
	const bool selectAllOfType = renderContext.selectedMolecule
		&& renderContext.selectedMolecule->name == molecule->name
		&& std::ranges::find(renderContext.selectedMolecule->atomIds, atomId) != renderContext.selectedMolecule->atomIds.end()
		&& lastSelectedAtomId == atomId
		&& renderContext.selectedMolecule->number > 0;
	if (selectAllOfType) {
		Rendering::MoleculeInfo selection{ molecule->name, 0, molecule->typeCount };
		for (const auto& candidate : molecules)
			if (candidate.name == molecule->name)
				selection.atomIds.insert(selection.atomIds.end(), candidate.atomIds.begin(), candidate.atomIds.end());
		renderContext.selectedMolecule = std::move(selection);
	}
	else {
		renderContext.selectedMolecule = *molecule;
	}
	_UpdateSelection(renderContext, std::set<int>{ renderContext.selectedMolecule->atomIds.begin(),
		renderContext.selectedMolecule->atomIds.end() });
	lastSelectedAtomId = atomId;
}

void Display::OnMouseButton(int button, int action, int mods) {
	if (!activeRenderContext)
		return;
	auto& renderContext = *activeRenderContext;
	if (ImGui::GetCurrentContext() && ImGui::GetIO().WantCaptureMouse) {
		if (action == GLFW_RELEASE) {
			isDragging = false;
			if (renderContext.activeGizmo) {
				renderContext.activeGizmo->activeAxis.reset();
				renderContext.activeGizmo->pullForce.reset();
				renderContext.activeGizmo->rotateForce.reset();
			}
		}
		return;
	}
	if (action == GLFW_PRESS)
		renderContext.revolveCamera = false;

	glm::ivec2 pixel{ static_cast<int>(mousePos.x),        static_cast<int>(mousePos.y) };
	const int objectId = GetObjectIdAtPixel(pixel);

    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        mousePosAtRightBtnDown = std::nullopt;
        if (action == GLFW_PRESS) {
            glfwGetCursorPos(window, &mousePos.x, &mousePos.y);
            mousePosAtBtnDown = mousePos;
            timeAtBtnDown = std::chrono::steady_clock::now();
            isDragging = true;

			if (gizmoEnabled && renderContext.activeGizmo) {
				renderContext.activeGizmo->SetActiveAxis(objectId);
				renderContext.activeGizmo->BeginDragging(mousePos, *renderContext.camera, windowSize);
            }

        }
        else if (action == GLFW_RELEASE) {
			if (renderContext.activeGizmo) {
				renderContext.activeGizmo->activeAxis.reset();
				renderContext.activeGizmo->pullForce.reset();
				stopMovingLiveeditCmd.store(true);
            }

            isDragging = false;

            glm::dvec2 mousePos{};
            glfwGetCursorPos(window, &mousePos.x, &mousePos.y);
            auto durationMs = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - timeAtBtnDown).count();
            bool isClick = glm::distance(mousePos, mousePosAtBtnDown) < 5. && durationMs < 200;

            if (isClick && drawAtomsFromCpuShader) {
                if (allowUserInputs) {
                    bool objectIdIsAtomId = ElementIdIsAtomid(objectId);
                    int id = objectIdIsAtomId ? objectId : -1;
                    std::lock_guard<std::mutex> lock(liveEditCommandsQueueMutex);
                    liveEditCommandsQueue.push_back(LiveEdit::AtomSelected{ id });
                    HandleGizmo(objectId);
                }
                else if (ElementIdIsAtomid(objectId))
                    SelectMolecule(objectId);
                else
                    SelectMolecule(-1);
            }
        }
    }
    else if (button == GLFW_MOUSE_BUTTON_RIGHT) {
        if (action == GLFW_PRESS) {
            mousePosAtRightBtnDown.emplace();
            glfwGetCursorPos(window, &mousePosAtRightBtnDown->x, &mousePosAtRightBtnDown->y);            
		}
    }
}

void Display::OnMouseScroll(double xoffset, double yoffset) {
	if (ImGui::GetCurrentContext() && ImGui::GetIO().WantCaptureMouse)
		return;
	if (activeRenderContext)
		activeRenderContext->camera->Update(0, 0, yoffset * 0.1f);
}
// -------------------------------------------------------------------------------------------------- //


void Display::ConsumeInputs(bool& shouldRecolorAtoms) {
	if (!activeRenderContext)
		return;
	auto& renderContext = *activeRenderContext;

	if (gizmoEnabled && renderContext.activeGizmo) {
		renderContext.activeGizmo->UpdateDraggingForce(mousePos, *renderContext.camera, windowSize);
    }


    while (!overlay->submittedCommands.empty()) {
        std::visit([&](auto&& cmd) {
            using T = std::decay_t<decltype(cmd)>;
            if constexpr (std::is_same_v<T, Overlay::SubmittedCmd>) {
                std::optional<LiveEdit::Command> liveEditCmd = LiveEdit::ParseCommand(cmd.cmd);
                if (liveEditCmd.has_value()) {
                    std::lock_guard<std::mutex> lock(liveEditCommandsQueueMutex);
                    liveEditCommandsQueue.push_back(*liveEditCmd);
                }
            }
            else if constexpr (std::is_same_v<T, ColoringMethod>) {
				renderContext.renderSettings->coloringMethod = cmd;
                shouldRecolorAtoms = true;
            }
            else if constexpr (std::is_same_v<T, Overlay::SolventVisibility>) {
				renderContext.renderSettings->showSolvents = cmd.visible;
                shouldRecolorAtoms |= true;
            }
			else if constexpr (std::is_same_v<T, Overlay::ResetCamera>) {
				renderContext.camera->Reset();
				renderContext.revolveCamera = false;
			}
			else if constexpr (std::is_same_v<T, Overlay::RevolveCamera>) {
				renderContext.revolveCamera = true;
				renderContext.lastRevolveTime = std::chrono::high_resolution_clock::now();
			}
			else if constexpr (std::is_same_v<T, Overlay::SelectSimulation>) {
				activeSimulationId = cmd.simulationId;
				isDragging = false;
				mousePosAtRightBtnDown.reset();
			}
            else {
                int a = 0;
            }
			}, overlay->submittedCommands.front());
		overlay->submittedCommands.pop_front();
    }
}
