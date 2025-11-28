#include <GL/glew.h>
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "backends/imgui_impl_opengl3.h"

#include "Display.h"
#include "filesystem"


static void RightAlignedField(const std::string& label,
    const std::string& value,
    const std::string& unit,
    const std::string& maxPattern)
{
    // Label
    ImGui::Text("%s", label.c_str());
    ImGui::SameLine();

    // Width calculations
    const ImVec2 maxWidth = ImGui::CalcTextSize(maxPattern.c_str());
    const ImVec2 valWidth = ImGui::CalcTextSize(value.c_str());

    float pad = maxWidth.x - valWidth.x;
    if (pad < 0.f) pad = 0.f;

    // pad to right-align
    ImGui::SetCursorPosX(ImGui::GetCursorPosX() + pad);

    // Value + optional unit
    if (unit.empty()) {
        ImGui::Text("%s ", value.c_str());
    }
    else {
        ImGui::Text("%s %s ", value.c_str(), unit.c_str());
    }

    ImGui::SameLine();   // keep next field on same row
}









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

void DrawTopBar(const SimStatus& status) {
    const float barHeight = 36.0f;
    ImGuiIO& io = ImGui::GetIO();

    ImGui::SetNextWindowPos(ImVec2(0, 0));
    ImGui::SetNextWindowSize(ImVec2(io.DisplaySize.x, barHeight));

    ImGuiWindowFlags flags =
        ImGuiWindowFlags_NoTitleBar
        | ImGuiWindowFlags_NoResize
        | ImGuiWindowFlags_NoMove
        | ImGuiWindowFlags_NoScrollbar
        | ImGuiWindowFlags_NoSavedSettings
        | ImGuiWindowFlags_NoNav
        | ImGuiWindowFlags_NoBringToFrontOnFocus;

    ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.0f);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.0f);
    ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0.10f, 0.10f, 0.12f, 0.95f));

    ImGui::Begin("###TopStatusBar", nullptr, flags);

    ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(14, 0));
    ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(0, 0));

    RightAlignedField(
        "Step",
        std::to_string(status.step),
        "",
        "999999999"   // width reference
    );

    // Temperature
    if (status.temperature) {
        RightAlignedField(
            "Temp",
            std::format("{:.2f}", *status.temperature),
            "[K]",
            "9999.99"
        );
    }

    // Max force
    if (status.maxForce) {
        RightAlignedField(
            "MaxF",
            std::format("{:.2f}", *status.maxForce),
            "[kJ/mol/nm]",
            "99999999.99"
        );
    }

    // Performance
    RightAlignedField(
        "Performance",
        std::format("{:.3f}", status.avgStepTime),
        "[ms/step]",
        "999.999"
    );
    if (status.simulationPerformance) {
        RightAlignedField(
            "",
            std::format("{:.2f}", *status.simulationPerformance),
            "[ns/day]",
            "999.99"
        );
    }


    ImGui::PopStyleVar(2);
    ImGui::End();

    ImGui::PopStyleVar(2); // rounding + border size
    ImGui::PopStyleColor(); // bg
}

void DrawBottomBar(RenderSettings& renderSettings) {

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

void Overlay::Draw(RenderSettings& renderSettings, const SimStatus& simstatus) {
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();

	DrawTopBar(simstatus);
    DrawBottomBar(renderSettings);
}


void Overlay::Render() {
    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}