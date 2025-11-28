#include <GL/glew.h>
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "backends/imgui_impl_opengl3.h"

#include "Display.h"
#include "filesystem"

Overlay::Overlay(GLFWwindow* window, const std::filesystem::path& limaDir) {
    // Setup imgui
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
    io.Fonts->AddFontFromFileTTF(
		(limaDir / "resources" / "ui" / "Roboto-Medium.ttf").string().c_str(),
        22.0f
    );

    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 430");

}
Overlay::~Overlay() {
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
}

void Overlay::Draw(RenderSettings& renderSettings) {
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();

    float barHeight = 50.0f;
    ImVec2 winSize = ImGui::GetIO().DisplaySize;

    ImGui::SetNextWindowPos(ImVec2(0, winSize.y - barHeight));
    ImGui::SetNextWindowSize(ImVec2(winSize.x, barHeight));

    ImGuiWindowFlags flags =
        ImGuiWindowFlags_NoTitleBar
        | ImGuiWindowFlags_NoResize
        | ImGuiWindowFlags_NoMove
        | ImGuiWindowFlags_NoCollapse
        | ImGuiWindowFlags_NoScrollbar
        | ImGuiWindowFlags_NoSavedSettings
        | ImGuiWindowFlags_NoBringToFrontOnFocus
        | ImGuiWindowFlags_NoNav;

    ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.0f);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.0f);
    ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0.12f, 0.12f, 0.12f, 0.95f));

    ImGui::Begin("BottomBar", nullptr, flags);

    ImGui::Checkbox("Show solvents", &renderSettings.showSolvents);
    ImGui::SameLine();
    //ImGui::Checkbox("Render facets", &renderSettings.renderFacets);
    //ImGui::SameLine();
    //ImGui::Checkbox("Render normals", &renderSettings.renderFacetNormals);

    ImGui::End();

    ImGui::PopStyleVar(2);
    ImGui::PopStyleColor();
}


void Overlay::Render() {
    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}