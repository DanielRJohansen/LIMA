#include <GL/glew.h>
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "backends/imgui_impl_opengl3.h"

#include "Display.h"
#include "filesystem"

#include <deque>
#include <string>
#include <cstdio>
#include <format>


void RightAlignedField(const std::string& label,
    const std::string& value,
    const std::string& unit,
    const std::string& maxPattern)
{
    ImGui::Text("%s", label.c_str());
    ImGui::SameLine();

    const ImVec2 maxWidth = ImGui::CalcTextSize(maxPattern.c_str());
    const ImVec2 valWidth = ImGui::CalcTextSize(value.c_str());

    float pad = maxWidth.x - valWidth.x;
    if (pad < 0.f) pad = 0.f;

    ImGui::SetCursorPosX(ImGui::GetCursorPosX() + pad);

    if (unit.empty()) {
        ImGui::Text("%s ", value.c_str());
    }
    else {
        ImGui::Text("%s %s ", value.c_str(), unit.c_str());
    }

    ImGui::SameLine();
}

static std::deque<std::string>& ConsoleLines()
{
    static std::deque<std::string> lines;
    return lines;
}

char* ConsoleInputBuffer()
{
    static char buffer[512] = "";
    return buffer;
}

constexpr const char* ConsolePrompt()
{
    return "> ";
}

std::string SubmitConsoleInput()
{
    char* buffer = ConsoleInputBuffer();
    if (buffer[0] == '\0')
        return "";

    auto& lines = ConsoleLines();
    lines.emplace_back(std::format("{}{}", ConsolePrompt(), buffer));
    if (lines.size() > 2)
        lines.pop_front();
	std::string inputText(buffer);


    //std::printf("[OverlayConsole] %s\n", buffer);
    buffer[0] = '\0';
	return inputText;
}

int TerminalInputCallback(ImGuiInputTextCallbackData* data)
{
    if (data->EventFlag == ImGuiInputTextFlags_CallbackHistory) {
        return 0;
    }
    return 0;
}

Overlay::Overlay(GLFWwindow* window, const std::filesystem::path& limaDir) {
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
    io.Fonts->AddFontFromFileTTF(
        (limaDir / "resources" / "ui" / "Roboto-Medium.ttf").string().c_str(),
        22.0f
    );
    io.IniFilename = nullptr; // disable imgui.ini creation

    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 430");
}

Overlay::~Overlay() {
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
}

void DrawTopBar(const SimStatus& status, int fps) {
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

    RightAlignedField("Step", std::to_string(status.step), "", "999999999");

    if (status.temperature) {
        RightAlignedField(
            "Temp",
            std::format("{:.2f}", *status.temperature),
            "[K]",
            "9999.99"
        );
    }

    if (status.maxForce) {
        RightAlignedField(
            "MaxF",
            std::format("{:.2f}", *status.maxForce),
            "[kJ/mol/nm]",
            "99999999.99"
        );
    }

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

#ifdef _DEBUG
    RightAlignedField("FPS", std::to_string(fps), "", "9999");
#endif

    ImGui::PopStyleVar(2);
    ImGui::End();

    ImGui::PopStyleVar(2);
    ImGui::PopStyleColor();
}

void Overlay::HandleConsole()
{
    constexpr float bottomBarHeight = 50.0f;
    constexpr float consoleHeight = 100.0f;

    const ImVec2 winSize = ImGui::GetIO().DisplaySize;

    ImGui::SetNextWindowPos(ImVec2(0, winSize.y - bottomBarHeight - consoleHeight));
    ImGui::SetNextWindowSize(ImVec2(winSize.x, consoleHeight));

    ImGuiWindowFlags flags =
        ImGuiWindowFlags_NoTitleBar
        | ImGuiWindowFlags_NoResize
        | ImGuiWindowFlags_NoMove
        | ImGuiWindowFlags_NoCollapse
        | ImGuiWindowFlags_NoSavedSettings
        | ImGuiWindowFlags_NoBringToFrontOnFocus
        | ImGuiWindowFlags_NoScrollbar;

    ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.0f);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.0f);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(8.0f, 8.0f));
    ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0.05f, 0.05f, 0.06f, 0.97f));
    ImGui::PushStyleColor(ImGuiCol_FrameBg, ImVec4(0.05f, 0.05f, 0.06f, 0.0f));
    ImGui::PushStyleColor(ImGuiCol_Border, ImVec4(0.0f, 0.0f, 0.0f, 0.0f));

    ImGui::Begin("OverlayConsole", nullptr, flags);

    auto& lines = ConsoleLines();
    char* inputBuffer = ConsoleInputBuffer();

    for (const std::string& line : lines) {
        ImGui::TextUnformatted(line.c_str());
    }

    ImGui::TextUnformatted(ConsolePrompt());
    ImGui::SameLine(0.0f, 0.0f);

    ImGui::PushItemWidth(-1.0f);
    const bool submitted = ImGui::InputText(
        "##TerminalInput",
        inputBuffer,
        512,
        ImGuiInputTextFlags_EnterReturnsTrue
        | ImGuiInputTextFlags_CallbackHistory,
        TerminalInputCallback
    );
    ImGui::PopItemWidth();

    if (ImGui::IsWindowAppearing()) {
        ImGui::SetKeyboardFocusHere(-1);
    }

    if (submitted) {
        std::string submittedCommand = SubmitConsoleInput();
		std::lock_guard<std::mutex> lock(consoleMutex);
		submittedCommands.push_back(submittedCommand);
        ImGui::SetKeyboardFocusHere(-1);
    }

    ImGui::SetScrollHereY(1.0f);

    ImGui::End();

    ImGui::PopStyleColor(3);
    ImGui::PopStyleVar(3);
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

    ImGui::End();

    ImGui::PopStyleVar(2);
    ImGui::PopStyleColor();
}

void Overlay::Draw(RenderSettings& renderSettings, const SimStatus& simstatus, int fps) {
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();

    DrawTopBar(simstatus, fps);
    HandleConsole();
    DrawBottomBar(renderSettings);
    didDrawThisFrame = true;
}

void Overlay::Render() {
    if (!didDrawThisFrame)
        return;
    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
    didDrawThisFrame = false;
}