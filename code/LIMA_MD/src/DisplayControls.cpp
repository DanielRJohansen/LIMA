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












std::optional<int> GetHoveredGizmoAxisAtPixel(std::optional<TranslateGizmo>& gizmo, glm::ivec2 pixel) {
    if (!gizmo.has_value())
        return std::nullopt;

    // TODO: replace with your actual gizmo picking
    return std::nullopt;
}

void Display::OnMouseMove(double xpos, double ypos) {
    if (activeGizmo.has_value()) {
        activeGizmo->hoveredAxis = GetHoveredGizmoAxisAtPixel(activeGizmo, glm::ivec2{ (int)xpos, (int)ypos });
    }

    if (activeGizmo.has_value() && activeGizmo->isDragging) {
        mousePos.x = xpos;
        mousePos.y = ypos;
        return;
    }

    if (isDragging) {
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
    if (atomId != -1)
        renderAtomsHost[atomId].HighLight(true);
    prevAtomId = atomId;
    renderAtoms.SetData(renderAtomsHost);
}

void Display::HandleGizmo(int atomId) {
    if (atomId == -1) {
        activeGizmo.reset();
        return;
    }

    if (!activeGizmo.has_value()) {
        activeGizmo = TranslateGizmo{};
    }

    activeGizmo->position = glm::vec3{ renderAtomsTemp[atomId].position.x, renderAtomsTemp[atomId].position.y, renderAtomsTemp[atomId].position.z };
}



void Display::OnMouseButton(int button, int action, int mods) {
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        if (action == GLFW_PRESS) {
            glfwGetCursorPos(window, &mousePos.x, &mousePos.y);
            mousePosAtBtnDown = mousePos;
            timeAtBtnDown = std::chrono::steady_clock::now();

            if (activeGizmo.has_value()) {
                activeGizmo->hoveredAxis = GetHoveredGizmoAxisAtPixel(activeGizmo, glm::ivec2{ (int)mousePos.x, (int)mousePos.y });

                if (activeGizmo->hoveredAxis.has_value()) {
                    activeGizmo->isDragging = true;
                    activeGizmo->activeAxis = activeGizmo->hoveredAxis;
                    activeGizmo->dragStartPosition = activeGizmo->position;
                    activeGizmo->dragStartMousePos = mousePos;
                    return;
                }
            }

            isDragging = true;
        }
        else if (action == GLFW_RELEASE) {
            if (activeGizmo.has_value() && activeGizmo->isDragging) {
                activeGizmo->isDragging = false;
                activeGizmo->activeAxis.reset();
                activeGizmo->hoveredAxis = GetHoveredGizmoAxisAtPixel(activeGizmo, glm::ivec2{ (int)mousePos.x, (int)mousePos.y });
                return;
            }

            isDragging = false;

            glm::dvec2 mousePos{};
            glfwGetCursorPos(window, &mousePos.x, &mousePos.y);
            auto durationMs = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - timeAtBtnDown).count();
            bool isClick = glm::distance(mousePos, mousePosAtBtnDown) < 5. && durationMs < 200;

            if (isClick && drawAtomsFromCpuShader) {
                drawAtomsFromCpuShader->DrawPicking(*renderTargetControl);
                int atomId = renderTargetControl->ReadIdAtPixel(glm::ivec2{ (int)mousePos.x, (int)mousePos.y });
                HandleHighlightAtom(atomId, lastSelectedAtomId, drawAtomsFromCpuShader->renderAtomsBuffer);
                HandleGizmo(atomId);
            }
        }
    }
}

void Display::OnMouseScroll(double xoffset, double yoffset) {
    camera.Update(0, 0, yoffset * 0.1f);
}