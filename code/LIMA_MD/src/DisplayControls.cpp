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













void Display::OnMouseMove(double xpos, double ypos) {
    /*if (activeGizmo.has_value()) {
        activeGizmo->hoveredAxis = GetHoveredGizmoAxisAtPixel(activeGizmo, glm::ivec2{ (int)xpos, (int)ypos });
    }*/

    if (activeGizmo.has_value() && activeGizmo->activeAxis.has_value()) {
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