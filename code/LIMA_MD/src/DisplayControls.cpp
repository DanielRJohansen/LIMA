#include "Display.h"

#include <GLFW/glfw3.h>
#include <algorithm>




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
    if (activeGizmo.has_value() && activeGizmo->activeAxis.has_value()) {
        activeGizmo->UpdateDraggingForce(glm::vec2(xpos, ypos), camera, windowSize);
    }
    else if (isDragging) {
        const float sensitivity = 0.001f;
        const float xOffset = static_cast<float>(xpos - mousePos.x) * sensitivity;
        const float yOffset = static_cast<float>(mousePos.y - ypos) * sensitivity;

        camera.Update(xOffset, -yOffset, 0);
    }

    mousePos.x = xpos;
    mousePos.y = ypos;
}

void Display::HandleGizmo(int objectId) {
    if (!allowUserInputs)
        return;

    if (objectId == -1) {
        activeGizmo.reset();
        stopMovingLiveeditCmd.store(true);
        return;
    }

    if (!activeGizmo.has_value()) {
        activeGizmo = TransformGizmo{};
    }
    if (objectId < renderAtomsHost.size()) {
        activeGizmo->position = glm::vec3{ renderAtomsHost[objectId].position.x, renderAtomsHost[objectId].position.y, renderAtomsHost[objectId].position.z };
        activeGizmo->idOfAtomAttachedTo = objectId;
    }
}

void Display::OnMouseButton(int button, int action, int mods) {

	glm::ivec2 pixel{ static_cast<int>(mousePos.x),        static_cast<int>(mousePos.y) };
	const int objectId = GetObjectIdAtPixel(pixel);

    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        mousePosAtRightBtnDown = std::nullopt;
        if (action == GLFW_PRESS) {
            glfwGetCursorPos(window, &mousePos.x, &mousePos.y);
            mousePosAtBtnDown = mousePos;
            timeAtBtnDown = std::chrono::steady_clock::now();
            isDragging = true;

            if (activeGizmo) {
				activeGizmo->SetActiveAxis(objectId);
                activeGizmo->BeginDragging(mousePos, camera, windowSize);
            }

        }
        else if (action == GLFW_RELEASE) {
            if (activeGizmo.has_value()) {
                activeGizmo->activeAxis.reset();
                activeGizmo->pullForce.reset();
				stopMovingLiveeditCmd.store(true);
            }

            isDragging = false;

            glm::dvec2 mousePos{};
            glfwGetCursorPos(window, &mousePos.x, &mousePos.y);
            auto durationMs = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - timeAtBtnDown).count();
            bool isClick = glm::distance(mousePos, mousePosAtBtnDown) < 5. && durationMs < 200;

            if (isClick && drawAtomsFromCpuShader) {
                {
                    bool objectIdIsAtomId = ElementIdIsAtomid(objectId);
                    int id = objectIdIsAtomId ? objectId : -1;
                    std::lock_guard<std::mutex> lock(liveEditCommandsQueueMutex);
                    liveEditCommandsQueue.push_back(LiveEdit::AtomSelected{ id });
                }

                HandleGizmo(objectId);
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
    camera.Update(0, 0, yoffset * 0.1f);
}
// -------------------------------------------------------------------------------------------------- //


void Display::ConsumeInputs(bool& shouldRecolorAtoms) {

    if (activeGizmo) {
        activeGizmo->UpdateDraggingForce(mousePos, camera, windowSize);
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
                rendersettings.coloringMethod = cmd;
                shouldRecolorAtoms |= true;
            }
            else if constexpr (std::is_same_v<T, Overlay::SolventVisibility>) {
                rendersettings.showSolvents = cmd.visible;
                shouldRecolorAtoms |= true;
            }
            else {
                int a = 0;
            }
			}, overlay->submittedCommands.front());
		overlay->submittedCommands.pop_front();
    }
}
