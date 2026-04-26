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

namespace
{
    constexpr ImVec4 kPanelBg = ImVec4(0.10f, 0.11f, 0.13f, 0.78f);
    constexpr ImVec4 kPanelBgStrong = ImVec4(0.12f, 0.13f, 0.16f, 0.88f);
    constexpr ImVec4 kPanelBorder = ImVec4(0.28f, 0.34f, 0.40f, 0.30f);
    constexpr ImVec4 kText = ImVec4(0.88f, 0.92f, 0.96f, 1.00f);
    constexpr ImVec4 kTextDim = ImVec4(0.60f, 0.67f, 0.74f, 1.00f);
    constexpr ImVec4 kAccent = ImVec4(0.38f, 0.63f, 0.92f, 1.00f);
    constexpr ImVec4 kWidget = ImVec4(0.18f, 0.20f, 0.24f, 0.95f);
    constexpr ImVec4 kWidgetHover = ImVec4(0.23f, 0.26f, 0.31f, 0.95f);
    constexpr ImVec4 kWidgetActive = ImVec4(0.28f, 0.32f, 0.38f, 0.95f);

    constexpr float kOuterMargin = 18.0f;
    constexpr float kPanelRounding = 8.0f;
    constexpr float kPanelBorderSize = 1.0f;
    constexpr float kTopBarHeight = 56.0f;
    constexpr float kBottomBarHeight = 62.0f;
    constexpr float kConsoleHeight = 128.0f;

    void PushOverlayTheme()
    {
        ImGuiStyle& style = ImGui::GetStyle();

        style.WindowRounding = kPanelRounding;
        style.ChildRounding = 14.0f;
        style.FrameRounding = 12.0f;
        style.PopupRounding = 12.0f;
        style.GrabRounding = 12.0f;
        style.ScrollbarRounding = 12.0f;
        style.TabRounding = 12.0f;

        style.WindowBorderSize = kPanelBorderSize;
        style.FrameBorderSize = 0.0f;
        style.PopupBorderSize = 0.0f;
        style.TabBorderSize = 0.0f;

        style.WindowPadding = ImVec2(16.0f, 12.0f);
        style.FramePadding = ImVec2(12.0f, 9.0f);
        style.ItemSpacing = ImVec2(14.0f, 10.0f);
        style.ItemInnerSpacing = ImVec2(8.0f, 6.0f);

        ImVec4* colors = style.Colors;
        colors[ImGuiCol_Text] = kText;
        colors[ImGuiCol_TextDisabled] = kTextDim;

        colors[ImGuiCol_WindowBg] = kPanelBg;
        colors[ImGuiCol_ChildBg] = ImVec4(0.00f, 0.00f, 0.00f, 0.00f);
        colors[ImGuiCol_PopupBg] = kPanelBgStrong;
        colors[ImGuiCol_Border] = kPanelBorder;
        colors[ImGuiCol_BorderShadow] = ImVec4(0, 0, 0, 0);

        colors[ImGuiCol_FrameBg] = kWidget;
        colors[ImGuiCol_FrameBgHovered] = kWidgetHover;
        colors[ImGuiCol_FrameBgActive] = kWidgetActive;

        colors[ImGuiCol_TitleBg] = kPanelBgStrong;
        colors[ImGuiCol_TitleBgActive] = kPanelBgStrong;
        colors[ImGuiCol_TitleBgCollapsed] = kPanelBg;

        colors[ImGuiCol_Button] = kWidget;
        colors[ImGuiCol_ButtonHovered] = kWidgetHover;
        colors[ImGuiCol_ButtonActive] = kWidgetActive;

        colors[ImGuiCol_Header] = kWidget;
        colors[ImGuiCol_HeaderHovered] = kWidgetHover;
        colors[ImGuiCol_HeaderActive] = kWidgetActive;

        colors[ImGuiCol_CheckMark] = kAccent;
        colors[ImGuiCol_SliderGrab] = kAccent;
        colors[ImGuiCol_SliderGrabActive] = ImVec4(0.88f, 0.79f, 0.67f, 1.00f);

        colors[ImGuiCol_ScrollbarBg] = ImVec4(0.10f, 0.09f, 0.08f, 0.35f);
        colors[ImGuiCol_ScrollbarGrab] = ImVec4(0.38f, 0.34f, 0.30f, 0.80f);
        colors[ImGuiCol_ScrollbarGrabHovered] = ImVec4(0.46f, 0.41f, 0.36f, 0.90f);
        colors[ImGuiCol_ScrollbarGrabActive] = ImVec4(0.54f, 0.48f, 0.42f, 1.00f);
    }

    void DrawPanelShadow(const ImVec2& min, const ImVec2& max, float rounding)
    {
        ImDrawList* drawList = ImGui::GetBackgroundDrawList();
        drawList->AddRectFilled(
            ImVec2(min.x + 0.0f, min.y + 8.0f),
            ImVec2(max.x + 0.0f, max.y + 8.0f),
            IM_COL32(0, 0, 0, 55),
            rounding
        );
        drawList->AddRectFilled(
            ImVec2(min.x + 0.0f, min.y + 16.0f),
            ImVec2(max.x + 0.0f, max.y + 16.0f),
            IM_COL32(0, 0, 0, 20),
            rounding
        );
    }

    void BeginFloatingPanel(const char* name, const ImVec2& pos, const ImVec2& size, ImGuiWindowFlags flags, bool strongBg = false)
    {
        DrawPanelShadow(pos, ImVec2(pos.x + size.x, pos.y + size.y), kPanelRounding);

        ImGui::SetNextWindowPos(pos);
        ImGui::SetNextWindowSize(size);

        ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, kPanelRounding);
        ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, kPanelBorderSize);
        ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(16.0f, 12.0f));
        ImGui::PushStyleColor(ImGuiCol_WindowBg, strongBg ? kPanelBgStrong : kPanelBg);
        ImGui::PushStyleColor(ImGuiCol_Border, kPanelBorder);

        ImGui::Begin(name, nullptr, flags);
    }

    void EndFloatingPanel()
    {
        ImGui::End();
        ImGui::PopStyleColor(2);
        ImGui::PopStyleVar(3);
    }

    void RightAlignedField(const std::string& label, const std::string& value, const std::string& unit, const std::string& maxPattern)
    {
        if (!label.empty()) {
            ImGui::PushStyleColor(ImGuiCol_Text, kTextDim);
            ImGui::TextUnformatted(label.c_str());
            ImGui::PopStyleColor();
            ImGui::SameLine();
        }

        const ImVec2 maxWidth = ImGui::CalcTextSize(maxPattern.c_str());
        const std::string combined = unit.empty() ? value : std::format("{} {}", value, unit);
        const ImVec2 valueWidth = ImGui::CalcTextSize(combined.c_str());

        float pad = maxWidth.x - valueWidth.x;
        if (pad > 0.0f)
            ImGui::SetCursorPosX(ImGui::GetCursorPosX() + pad);

        if (unit.empty()) {
            ImGui::TextUnformatted(value.c_str());
        }
        else {
            ImGui::Text("%s %s", value.c_str(), unit.c_str());
        }

        ImGui::SameLine();
    }






    const char* ColoringMethodName(ColoringMethod coloringMethod)
    {
        switch (coloringMethod) {
        case ColoringMethod::Atomname: return "Atom name";
        case ColoringMethod::Charge: return "Charge";
        case ColoringMethod::GradientFromAtomid: return "Gradient from atom id";
        case ColoringMethod::PersistentClusterId: return "Gradient from compound id";
        default: return "Unknown";
        }
    }

    bool ColoringMethodMenuItem(
        std::deque<Overlay::Command>& submittedCommands,
        ColoringMethod currentMethod,
        ColoringMethod method
    ) {
        const bool isSelected = currentMethod == method;

        if (ImGui::MenuItem(ColoringMethodName(method), nullptr, isSelected)) {
            submittedCommands.push_back(method);
            return true;
        }

        return false;
    }

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
    buffer[0] = '\0';
    return inputText;
}

int TerminalInputCallback(ImGuiInputTextCallbackData* data)
{
    if (data->EventFlag == ImGuiInputTextFlags_CallbackHistory)
        return 0;
    return 0;
}

Overlay::Overlay(GLFWwindow* window, const std::filesystem::path& limaDir)
{
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();

    ImGuiIO& io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
    io.IniFilename = nullptr;

    // Better default choice than Roboto for this kind of UI.
    // Put Inter-Medium.ttf in resources/ui if you have it.
    if (std::filesystem::exists(limaDir / "resources" / "ui" / "Inter-Medium.ttf")) {
        io.Fonts->AddFontFromFileTTF(
            (limaDir / "resources" / "ui" / "Inter-Medium.ttf").string().c_str(),
            24.0f
        );
    }
    else {
        io.Fonts->AddFontFromFileTTF(
            (limaDir / "resources" / "ui" / "Roboto-Medium.ttf").string().c_str(),
            24.0f
        );
    }

    PushOverlayTheme();

    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 430");
}

Overlay::~Overlay()
{
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
}

void DrawSimstatusCard(const SimStatus& status, int fps)
{
    ImGuiIO& io = ImGui::GetIO();

    const bool hasSimulationStatus =
        status.step.has_value()
        || status.temperature.has_value()
        || status.maxForce.has_value()
        || status.expectedTimeToFinish.has_value();

    const bool hasEnginePerformance =
        status.avgStepTime.has_value()
        || status.simulationPerformance.has_value()
#ifdef _DEBUG
        || true
#endif
        ;

    if (!hasSimulationStatus && !hasEnginePerformance)
        return;

    auto DrawField = [](const char* label, const std::string& value, const char* unit = nullptr)
        {
            ImGui::TableNextRow();

            ImGui::TableSetColumnIndex(0);
            ImGui::PushStyleColor(ImGuiCol_Text, kTextDim);
            ImGui::TextUnformatted(label);
            ImGui::PopStyleColor();

            ImGui::TableSetColumnIndex(1);

            const float startX = ImGui::GetCursorPosX();
            const float colWidth = ImGui::GetColumnWidth();
            const float rightX = startX + colWidth;

            if (unit && unit[0] != '\0') {
                constexpr float gap = 6.0f;

                const float unitWidth = ImGui::CalcTextSize(unit).x;
                const float valueWidth = ImGui::CalcTextSize(value.c_str()).x;

                const float unitX = rightX - unitWidth;
                const float valueX = unitX - gap - valueWidth;

                ImGui::SetCursorPosX(valueX);
                ImGui::TextUnformatted(value.c_str());

                ImGui::SameLine(0.0f, gap);
                ImGui::SetCursorPosX(unitX);
                ImGui::TextUnformatted(unit);
            }
            else {
                const float valueWidth = ImGui::CalcTextSize(value.c_str()).x;
                ImGui::SetCursorPosX(rightX - valueWidth);
                ImGui::TextUnformatted(value.c_str());
            }
        };
    int nRows = 0;
    if (hasSimulationStatus) {
        nRows += 1;
        if (status.step.has_value()) ++nRows;
        if (status.temperature.has_value()) ++nRows;
        if (status.maxForce.has_value()) ++nRows;
        if (status.expectedTimeToFinish.has_value()) ++nRows;
    }
    if (hasEnginePerformance) {
        if (nRows > 0)
            nRows += 1;
        nRows += 1;
        if (status.avgStepTime.has_value()) ++nRows;
        if (status.simulationPerformance.has_value()) ++nRows;
#ifdef _DEBUG
        ++nRows;
#endif
    }

    const float cardWidth = 380.0f;
    const float lineHeight = ImGui::GetTextLineHeight();
    const float verticalPadding = 10.0f;
    const float rowSpacing = 8.0f;
    const float titleSpacing = 10.0f;
    const float sectionSpacing = 12.0f;
    const float cardHeight =
        verticalPadding * 2.0f
        + nRows * lineHeight
        + (nRows - 1) * rowSpacing
        + titleSpacing;

    const ImVec2 pos(kOuterMargin, kOuterMargin);
    const ImVec2 size(cardWidth, cardHeight);

    ImGuiWindowFlags flags =
        ImGuiWindowFlags_NoTitleBar
        | ImGuiWindowFlags_NoResize
        | ImGuiWindowFlags_NoMove
        | ImGuiWindowFlags_NoScrollbar
        | ImGuiWindowFlags_NoSavedSettings
        | ImGuiWindowFlags_NoNav
        | ImGuiWindowFlags_NoBringToFrontOnFocus;

    BeginFloatingPanel("###TopStatusCard", pos, size, flags, true);

    ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(10.0f, rowSpacing));

    if (ImGui::BeginTable("##TopStatusTable", 2, ImGuiTableFlags_SizingStretchProp))
    {
        ImGui::TableSetupColumn("Label", ImGuiTableColumnFlags_WidthStretch);
        ImGui::TableSetupColumn("Value", ImGuiTableColumnFlags_WidthStretch);

        if (hasSimulationStatus) {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::PushStyleColor(ImGuiCol_Text, kAccent);
            ImGui::TextUnformatted("Simulation");
            ImGui::PopStyleColor();
            ImGui::TableSetColumnIndex(1);
            ImGui::Dummy(ImVec2(0.0f, 0.0f));

            if (status.step.has_value())
                DrawField("Step", std::to_string(*status.step));

            if (status.temperature.has_value())
                DrawField("Temperature", std::format("{:.2f}", *status.temperature), "[K]");

            if (status.maxForce.has_value())
                DrawField("Max force", std::format("{:.2e}", *status.maxForce), "[kJ/mol/nm]");

            if (status.expectedTimeToFinish.has_value())
                DrawField("Remaining time", StringUtils::FormatTime(*status.expectedTimeToFinish, 3, 2));
        }

        if (hasEnginePerformance) {
            if (hasSimulationStatus) {
                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0);
                ImGui::Dummy(ImVec2(0.0f, sectionSpacing));
                ImGui::TableSetColumnIndex(1);
                ImGui::Dummy(ImVec2(0.0f, sectionSpacing));
            }

            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::PushStyleColor(ImGuiCol_Text, kAccent);
            ImGui::TextUnformatted("Performance");
            ImGui::PopStyleColor();
            ImGui::TableSetColumnIndex(1);
            ImGui::Dummy(ImVec2(0.0f, 0.0f));

            if (status.avgStepTime.has_value())
                DrawField("Step time", std::format("{:.3f}", *status.avgStepTime), "[ms]");

            if (status.simulationPerformance.has_value())
                DrawField("Simulation", std::format("{:.2f}", *status.simulationPerformance), "[ns/day]");

#ifdef _DEBUG
            DrawField("FPS", std::to_string(fps));
#endif
        }

        ImGui::EndTable();
    }

    ImGui::PopStyleVar();
    EndFloatingPanel();
}

void Overlay::HandleConsole()
{
    const ImVec2 winSize = ImGui::GetIO().DisplaySize;

    const ImVec2 pos(
        kOuterMargin,
        winSize.y - kOuterMargin - kBottomBarHeight - 10.0f - kConsoleHeight
    );
    const ImVec2 size(
        winSize.x - 2.0f * kOuterMargin,
        kConsoleHeight
    );

    ImGuiWindowFlags flags =
        ImGuiWindowFlags_NoTitleBar
        | ImGuiWindowFlags_NoResize
        | ImGuiWindowFlags_NoMove
        | ImGuiWindowFlags_NoCollapse
        | ImGuiWindowFlags_NoSavedSettings
        | ImGuiWindowFlags_NoBringToFrontOnFocus
        | ImGuiWindowFlags_NoScrollbar;

    BeginFloatingPanel("OverlayConsole", pos, size, flags, true);

    auto& lines = ConsoleLines();
    char* inputBuffer = ConsoleInputBuffer();

    ImGui::PushStyleColor(ImGuiCol_Text, kTextDim);
    for (const std::string& line : lines)
        ImGui::TextUnformatted(line.c_str());
    ImGui::PopStyleColor();

    ImGui::Spacing();

    ImGui::PushStyleColor(ImGuiCol_Text, kAccent);
    ImGui::TextUnformatted(ConsolePrompt());
    ImGui::PopStyleColor();
    ImGui::SameLine(0.0f, 8.0f);

    ImGui::PushItemWidth(-1.0f);
    const bool submitted = ImGui::InputText(
        "##TerminalInput",
        inputBuffer,
        512,
        ImGuiInputTextFlags_EnterReturnsTrue | ImGuiInputTextFlags_CallbackHistory,
        TerminalInputCallback
    );
    ImGui::PopItemWidth();

    if (ImGui::IsWindowAppearing())
        ImGui::SetKeyboardFocusHere(-1);

    if (submitted) {
        SubmittedCmd submittedCommand{ SubmitConsoleInput() };
        submittedCommands.push_back(submittedCommand);
        ImGui::SetKeyboardFocusHere(-1);
    }

    EndFloatingPanel();
}

void Overlay::HandleContextMenu(RenderSettings& renderSettings, std::optional<glm::dvec2> rightClickedPos)
{
    if (rightClickedPos.has_value()) {
        ImGui::SetNextWindowPos(ImVec2(
            static_cast<float>(rightClickedPos->x),
            static_cast<float>(rightClickedPos->y)
        ));
        ImGui::OpenPopup("OverlayContextMenu");
    }

    if (ImGui::BeginPopup("OverlayContextMenu")) {
        ImGui::TextUnformatted("Coloring method");
        ImGui::Separator();

        ColoringMethodMenuItem(submittedCommands, renderSettings.coloringMethod, ColoringMethod::Atomname);
        ColoringMethodMenuItem(submittedCommands, renderSettings.coloringMethod, ColoringMethod::Charge);
        ColoringMethodMenuItem(submittedCommands, renderSettings.coloringMethod, ColoringMethod::GradientFromAtomid);
        ColoringMethodMenuItem(submittedCommands, renderSettings.coloringMethod, ColoringMethod::PersistentClusterId);

        ImGui::EndPopup();
    }
}

void DrawBottomBar(RenderSettings& renderSettings)
{
    ImVec2 winSize = ImGui::GetIO().DisplaySize;

    const ImVec2 pos(
        kOuterMargin,
        winSize.y - kOuterMargin - kBottomBarHeight
    );
    const ImVec2 size(
        winSize.x - 2.0f * kOuterMargin,
        kBottomBarHeight
    );

    ImGuiWindowFlags flags =
        ImGuiWindowFlags_NoTitleBar
        | ImGuiWindowFlags_NoResize
        | ImGuiWindowFlags_NoMove
        | ImGuiWindowFlags_NoCollapse
        | ImGuiWindowFlags_NoScrollbar
        | ImGuiWindowFlags_NoSavedSettings
        | ImGuiWindowFlags_NoBringToFrontOnFocus
        | ImGuiWindowFlags_NoNav;

    BeginFloatingPanel("BottomBar", pos, size, flags);

    float widgetHeight = ImGui::GetFrameHeight();
    float offset = (kBottomBarHeight - widgetHeight) * 0.5f;

    ImGui::SetCursorPosY(offset);
    ImGui::Checkbox("Show solvents", &renderSettings.showSolvents);

    EndFloatingPanel();
}

void Overlay::Draw(RenderSettings& renderSettings, const SimStatus& simstatus, int fps, std::optional<glm::dvec2> rightClickedPos)
{
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();

    DrawSimstatusCard(simstatus, fps);
    if (enableConsole)
        HandleConsole();

    DrawBottomBar(renderSettings);
    HandleContextMenu(renderSettings, rightClickedPos);

    didDrawThisFrame = true;
}

void Overlay::Render()
{
    if (!didDrawThisFrame)
        return;

    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
    didDrawThisFrame = false;
}