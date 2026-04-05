#include <GL/glew.h>
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "backends/imgui_impl_opengl3.h"


#include "Display.h"
#include "Shaders.h"    
#include "TimeIt.h"
#include "MDFiles.h"




#include <GLFW/glfw3.h>
#include <algorithm>




static glm::vec3 GetAxisDirection(int axis)
{
    switch (axis) {
    case 0: return glm::vec3(1.f, 0.f, 0.f);
    case 1: return glm::vec3(0.f, 1.f, 0.f);
    case 2: return glm::vec3(0.f, 0.f, 1.f);
    default: return glm::vec3(0.f);
    }
}

static glm::vec2 WorldToScreen(
    const glm::vec3& worldPos,
    const glm::mat4& view,
    const glm::mat4& proj,
    int viewportWidth,
    int viewportHeight)
{
    const glm::vec4 clip = proj * view * glm::vec4(worldPos, 1.f);
    if (clip.w == 0.f)
        return glm::vec2(0.f);

    const glm::vec3 ndc = glm::vec3(clip) / clip.w;

    return glm::vec2{
        (ndc.x * 0.5f + 0.5f) * static_cast<float>(viewportWidth),
        (1.f - (ndc.y * 0.5f + 0.5f)) * static_cast<float>(viewportHeight)
    };
}


void TranslateGizmo::SetActiveAxis(int selectedObjectId) {
    if (selectedObjectId == arrowX.uniqueId)
        activeAxis = 0;
    else if (selectedObjectId == arrowY.uniqueId)
        activeAxis = 1;
    else if (selectedObjectId == arrowZ.uniqueId)
        activeAxis = 2;
    else
        activeAxis = std::nullopt;
}
void TranslateGizmo::UpdateDraggingForce(glm::vec2 mousePos, const Camera& camera, glm::vec2 windowSize) {
    const glm::vec2 mouseDelta = mousePos - dragStartMousePos;

    // Replace these with your actual matrices/getters.
    const glm::mat4 view = camera.View();
    const glm::mat4 proj = camera.Projection();

    int viewport[4]{};
    glGetIntegerv(GL_VIEWPORT, viewport);
    const int viewportWidth = viewport[2];
    const int viewportHeight = viewport[3];

    const glm::vec3 axisDirWorld = GetAxisDirection(*activeAxis);
    const glm::vec2 gizmoScreenPos = WorldToScreen(position, view, proj, viewportWidth, viewportHeight);
    const glm::vec2 gizmoAxisScreenPos = WorldToScreen(position + axisDirWorld, view, proj, viewportWidth, viewportHeight);

    glm::vec2 axisDirScreen = gizmoAxisScreenPos - gizmoScreenPos;
    const float axisLen = glm::length(axisDirScreen);

    if (axisLen > 1e-5f) {
        axisDirScreen /= axisLen;

        const float signedPixels = glm::dot(mouseDelta, axisDirScreen);

        // 100 screen pixels = full pull.
        float pixelsForFullPull = (float)std::max(windowSize.x, windowSize.y);
        const float pull = std::clamp(signedPixels / pixelsForFullPull, -1.f, 1.f);

        glm::vec3 force{};
        force[activeAxis.value()] = pull;
		pullForce = force;
    }
}


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

void HandleHighlightAtom(int atomId, int& prevAtomId, SSBO& renderAtoms) {
    if (atomId == prevAtomId)
        return;

    auto renderAtomsHost = renderAtoms.GetData<RenderAtom>();
    if (prevAtomId != -1)
        renderAtomsHost[prevAtomId].HighLight(false);
    if (atomId != -1 && atomId < renderAtomsHost.size())
        renderAtomsHost[atomId].HighLight(true);

    prevAtomId = atomId < renderAtomsHost.size() ? atomId : -1;
    renderAtoms.SetData(renderAtomsHost);
}

void Display::HandleGizmo(int objectId) {
	//bool isGizmoElement = elementId == (int)UniqueRenderElementIds::gizmoArrowX || elementId == (int)UniqueRenderElementIds::gizmoArrowY || elementId == (int)UniqueRenderElementIds::gizmoArrowZ;
    if (objectId == -1) {
        activeGizmo.reset();
        return;
    }

    if (!activeGizmo.has_value()) {
        activeGizmo = TranslateGizmo{};
    }
    if (objectId < renderAtomsTemp.size()) {
        activeGizmo->position = glm::vec3{ renderAtomsTemp[objectId].position.x, renderAtomsTemp[objectId].position.y, renderAtomsTemp[objectId].position.z };
    }
}

void Display::OnMouseButton(int button, int action, int mods) {

	glm::ivec2 pixel{ static_cast<int>(mousePos.x),        static_cast<int>(mousePos.y) };
	const int objectId = GetObjectIdAtPixel(pixel);

    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        if (action == GLFW_PRESS) {
            glfwGetCursorPos(window, &mousePos.x, &mousePos.y);
            mousePosAtBtnDown = mousePos;
            timeAtBtnDown = std::chrono::steady_clock::now();
            isDragging = true;

            if (activeGizmo) {
				activeGizmo->SetActiveAxis(objectId);
                activeGizmo->dragStartPosition = activeGizmo->position;
                activeGizmo->dragStartMousePos = glm::vec2(mousePos);
            }

        }
        else if (action == GLFW_RELEASE) {
            if (activeGizmo.has_value()) {
                activeGizmo->activeAxis.reset();
            }

            isDragging = false;

            glm::dvec2 mousePos{};
            glfwGetCursorPos(window, &mousePos.x, &mousePos.y);
            auto durationMs = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - timeAtBtnDown).count();
            bool isClick = glm::distance(mousePos, mousePosAtBtnDown) < 5. && durationMs < 200;

            if (isClick && drawAtomsFromCpuShader) {
                HandleHighlightAtom(objectId, lastSelectedAtomId, drawAtomsFromCpuShader->renderAtomsBuffer);
                HandleGizmo(objectId);
            }
        }
    }
}

void Display::OnMouseScroll(double xoffset, double yoffset) {
    camera.Update(0, 0, yoffset * 0.1f);
}