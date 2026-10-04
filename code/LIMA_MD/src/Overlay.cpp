#include <GL/glew.h>
#include "Format.h"
#include <stb/stb_image.h>
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "backends/imgui_impl_opengl3.h"

#include "Display.h"
#include "DisplayInternal.h"
#include "Utilities.h"

#include <algorithm>
#include <array>
#include <filesystem>
#include <format>

namespace
{
    constexpr ImVec4 panelBg{ .055f, .063f, .078f, .96f };
    constexpr ImVec4 menuBg{ .035f, .041f, .051f, 1.f };
    constexpr ImVec4 tabBg{ .070f, .080f, .096f, 1.f };
    constexpr ImVec4 panelBorder{ .23f, .27f, .31f, .65f };
    constexpr ImVec4 textColor{ .89f, .92f, .94f, 1.f };
    constexpr ImVec4 mutedText{ .56f, .63f, .68f, 1.f };
    constexpr ImVec4 accent{ .48f, .91f, .76f, 1.f };
    constexpr ImVec4 selectedBg{ .10f, .22f, .20f, 1.f };
    constexpr float margin = 20.f;
    constexpr float consoleHeight = 164.f;
    constexpr ImGuiWindowFlags panelFlags = ImGuiWindowFlags_NoTitleBar
        | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove
        | ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoSavedSettings
        | ImGuiWindowFlags_NoBringToFrontOnFocus;

    enum class CardSide { Left, Right };

    struct CardLayout {
        float leftY;
        float rightY;
        float bottom;
        float left = 0.f;
        float width = 0.f;
        bool tiled = false;
        SimulationId simulationId = 0;
        float scale = 1.f;
    } cardLayout;

    struct CardLine {
        std::string text;
        std::string value;
        std::optional<std::string> unit;
        std::optional<float> progress;
    };

    struct CardSection {
        std::string title;
        std::vector<CardLine> lines;
    };

    void PushOverlayTheme()
    {
        auto& style = ImGui::GetStyle();
        style.WindowRounding = 10.f;
        style.ChildRounding = 6.f;
        style.FrameRounding = 5.f;
        style.PopupRounding = 8.f;
        style.GrabRounding = 4.f;
        style.ScrollbarRounding = 4.f;
        style.WindowBorderSize = 1.f;
        style.PopupBorderSize = 1.f;
        style.FrameBorderSize = 0.f;
        style.WindowPadding = ImVec2(18.f, 16.f);
        style.FramePadding = ImVec2(12.f, 7.f);
        style.ItemSpacing = ImVec2(10.f, 10.f);
        style.ItemInnerSpacing = ImVec2(8.f, 6.f);
        style.CellPadding = ImVec2(0.f, 6.f);
        style.ScrollbarSize = 8.f;

        auto* colors = style.Colors;
        colors[ImGuiCol_Text] = textColor;
        colors[ImGuiCol_TextDisabled] = mutedText;
        colors[ImGuiCol_WindowBg] = panelBg;
        colors[ImGuiCol_MenuBarBg] = menuBg;
        colors[ImGuiCol_PopupBg] = ImVec4(.065f, .078f, .09f, 1.f);
        colors[ImGuiCol_ChildBg] = ImVec4(0.f, 0.f, 0.f, 0.f);
        colors[ImGuiCol_Border] = panelBorder;
        colors[ImGuiCol_BorderShadow] = ImVec4(0.f, 0.f, 0.f, 0.f);
        colors[ImGuiCol_FrameBg] = ImVec4(.10f, .12f, .14f, 1.f);
        colors[ImGuiCol_FrameBgHovered] = ImVec4(.15f, .19f, .20f, 1.f);
        colors[ImGuiCol_FrameBgActive] = selectedBg;
        colors[ImGuiCol_Button] = ImVec4(.10f, .12f, .14f, 1.f);
        colors[ImGuiCol_ButtonHovered] = ImVec4(.16f, .22f, .23f, 1.f);
        colors[ImGuiCol_ButtonActive] = selectedBg;
        colors[ImGuiCol_Header] = selectedBg;
        colors[ImGuiCol_HeaderHovered] = ImVec4(.14f, .23f, .22f, 1.f);
        colors[ImGuiCol_HeaderActive] = selectedBg;
        colors[ImGuiCol_CheckMark] = accent;
        colors[ImGuiCol_SliderGrab] = accent;
        colors[ImGuiCol_SliderGrabActive] = accent;
        colors[ImGuiCol_Separator] = panelBorder;
        colors[ImGuiCol_ScrollbarBg] = ImVec4(0.f, 0.f, 0.f, 0.f);
        colors[ImGuiCol_ScrollbarGrab] = panelBorder;
        colors[ImGuiCol_ScrollbarGrabHovered] = mutedText;
        colors[ImGuiCol_ScrollbarGrabActive] = accent;
        colors[ImGuiCol_NavCursor] = accent;
        colors[ImGuiCol_TextSelectedBg] = selectedBg;
    }

    void BeginPanel(const char* name, ImVec2 pos, ImVec2 size, ImGuiWindowFlags flags = 0)
    {
        ImGui::SetNextWindowPos(pos);
        ImGui::SetNextWindowSize(size);
        ImGui::Begin(name, nullptr, panelFlags | flags);
    }

    void Metric(const char* label, const std::string& value, const char* unit);

    void DrawCard(const char* name, CardSide side, const std::vector<CardSection>& sections)
    {
        const float scale = cardLayout.scale;
        const float inset = margin * scale;
        const float width = std::max(1.f, std::min(380.f * scale, cardLayout.width - inset * 2.f));
        float& top = side == CardSide::Left || cardLayout.tiled ? cardLayout.leftY : cardLayout.rightY;
        if (cardLayout.bottom - top < 30.f * scale)
            return;
        const float x = cardLayout.left + (side == CardSide::Left || cardLayout.tiled ? inset : cardLayout.width - inset - width);
        ImGui::SetNextWindowSizeConstraints(ImVec2(width, 0.f), ImVec2(width, cardLayout.bottom - top));
        ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(14.f * scale, 10.f * scale));
        ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(8.f * scale, 4.f * scale));
        ImGui::PushStyleVar(ImGuiStyleVar_CellPadding, ImVec2(0.f, 2.f * scale));
        ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 10.f * scale);
        ImGui::PushStyleVar(ImGuiStyleVar_WindowMinSize, ImVec2(1.f, 1.f));
        ImGui::PushStyleVar(ImGuiStyleVar_ScrollbarSize, 8.f * scale);
        const auto windowName = Lima::Format("{}##{}", name, cardLayout.simulationId);
        BeginPanel(windowName.c_str(), ImVec2(x, top), ImVec2(width, 0.f), ImGuiWindowFlags_AlwaysAutoResize);
        for (size_t i = 0; i < sections.size(); ++i) {
            const auto& section = sections[i];
            if (i > 0) {
                ImGui::Spacing();
                ImGui::Separator();
                ImGui::Spacing();
            }
            ImGui::TextColored(accent, "%s", section.title.c_str());
            ImGui::PushID(static_cast<int>(i));
            bool tableOpen = false;
            for (const auto& line : section.lines) {
                if (line.progress) {
                    if (tableOpen) {
                        ImGui::EndTable();
                        tableOpen = false;
                    }
                    ImGui::PushStyleColor(ImGuiCol_PlotHistogram, accent);
                    ImGui::ProgressBar(*line.progress, ImVec2(-1.f, 5.f * scale), "");
                    ImGui::PopStyleColor();
                    continue;
                }
                if (!tableOpen) {
                    tableOpen = ImGui::BeginTable("##CardSection", 2, ImGuiTableFlags_SizingStretchProp);
                    if (tableOpen) {
                        ImGui::TableSetupColumn("Label", ImGuiTableColumnFlags_WidthStretch, .40f);
                        ImGui::TableSetupColumn("Value", ImGuiTableColumnFlags_WidthStretch, .60f);
                    }
                }
                if (tableOpen)
                    Metric(line.text.c_str(), line.value, line.unit ? line.unit->c_str() : "");
            }
            if (tableOpen)
                ImGui::EndTable();
            ImGui::PopID();
        }
        const float height = ImGui::GetWindowHeight();
        ImGui::End();
        ImGui::PopStyleVar(6);
        top += height + 12.f * scale;
    }

    void SectionLabel(const char* label)
    {
        ImGui::TextColored(mutedText, "%s", label);
        ImGui::Spacing();
    }

    void Tooltip(const char* text)
    {
        if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort))
            ImGui::SetTooltip("%s", text);
    }

    const char* ColoringMethodName(ColoringMethod method)
    {
        switch (method) {
        case ColoringMethod::Atomname: return "Atom name";
        case ColoringMethod::Charge: return "Charge";
        case ColoringMethod::GradientFromAtomid: return "Atom ID";
        case ColoringMethod::PersistentClusterId: return "Compound ID";
        case ColoringMethod::ForceMagnitude: return "Force magnitude";
        case ColoringMethod::NewCartoon: return "Backbone";
        default: return "Default";
        }
    }

    void ColoringMenu(const RenderSettings& settings, std::deque<Overlay::Command>& commands)
    {
        constexpr std::array methods{
            ColoringMethod::Atomname, ColoringMethod::Charge,
            ColoringMethod::GradientFromAtomid, ColoringMethod::PersistentClusterId,
            ColoringMethod::ForceMagnitude, ColoringMethod::NewCartoon
        };
        for (const auto method : methods) {
            if (method == ColoringMethod::ForceMagnitude && !settings.hasForceData)
                continue;
            if (method == ColoringMethod::NewCartoon && !settings.hasBackbone)
                continue;
            if (ImGui::MenuItem(ColoringMethodName(method), nullptr, settings.coloringMethod == method))
                commands.push_back(method);
        }
    }

    void CameraMenu(std::deque<Overlay::Command>& commands, SimulationId simulationId)
    {
        if (ImGui::MenuItem("Reset view"))
            commands.push_back(Overlay::ResetCamera{ simulationId });
        if (ImGui::MenuItem("Toggle orbit"))
            commands.push_back(Overlay::RevolveCamera{ simulationId });
    }

    void Metric(const char* label, const std::string& value, const char* unit = "")
    {
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0);
        ImGui::TextColored(mutedText, "%s", label);
        ImGui::TableSetColumnIndex(1);
        const float unitWidth = *unit ? ImGui::CalcTextSize(unit).x + 7.f * cardLayout.scale : 0.f;
        const float width = ImGui::CalcTextSize(value.c_str()).x + unitWidth;
        ImGui::SetCursorPosX(ImGui::GetCursorPosX() + std::max(0.f, ImGui::GetContentRegionAvail().x - width));
        ImGui::TextUnformatted(value.c_str());
        if (*unit) {
            ImGui::SameLine(0.f, 7.f * cardLayout.scale);
            ImGui::TextColored(mutedText, "%s", unit);
        }
    }

    void DrawTelemetry(const SimStatus& status)
    {
        if (!status.step && !status.temperature && !status.maxForce
            && !status.expectedTimeToFinish && !status.avgStepTime && !status.simulationPerformance)
            return;

        const bool hasSimulation = status.step || status.temperature || status.maxForce || status.expectedTimeToFinish;
        std::vector<CardSection> sections;
        if (hasSimulation) {
            auto& simulation = sections.emplace_back();
            simulation.title = "Simulation";
            if (status.step)
                simulation.lines.push_back({ "Step", Lima::Format("{}", *status.step) });
            if (status.progress)
                simulation.lines.push_back({ {}, {}, {}, std::clamp(*status.progress, 0.f, 1.f) });
            if (status.temperature)
                simulation.lines.push_back({ "Temperature", Lima::Format("{:.2f}", *status.temperature), "K" });
            if (status.maxForce)
                simulation.lines.push_back({ "Max force", Lima::Format("{:.2e}", *status.maxForce), "kJ/mol/nm" });
            if (status.expectedTimeToFinish)
                simulation.lines.push_back({ "Remaining", StringUtils::FormatTime(*status.expectedTimeToFinish, 3, 2) });
        }
        if (status.avgStepTime || status.simulationPerformance) {
            auto& performance = sections.emplace_back();
            performance.title = "Performance";
            if (status.avgStepTime)
                performance.lines.push_back({ "Step time", Lima::Format("{:.3f}", *status.avgStepTime), "ms" });
            if (status.simulationPerformance)
                performance.lines.push_back({ "Throughput", Lima::Format("{:.2f}", *status.simulationPerformance), "ns/day" });
        }
        DrawCard("##Telemetry", CardSide::Left, sections);
    }

    void DrawMoleculeInfo(const std::optional<Rendering::MoleculeInfo>& molecule)
    {
        if (!molecule)
            return;
		const std::string name = molecule->number > 0
			? Lima::Format("{} (#{})", molecule->name, molecule->number)
			: molecule->name;
		CardSection section{ "Molecule info", { { "Name", name } } };
		if (molecule->number == 0)
			section.lines.push_back({ "Molecules selected", Lima::Format("{}", molecule->typeCount) });
		section.lines.push_back({ "Atoms", Lima::Format("{}", molecule->atomIds.size()) });
        DrawCard("##MoleculeInfo", CardSide::Right, {
			std::move(section)
        });
    }

    float DrawMenuBar(RenderSettings& settings, std::deque<Overlay::Command>& commands, int fps,
        unsigned int logoTexture, bool tiled, bool enableConsole, SimulationId simulationId)
    {
        float height = 0.f;
        ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(12.f, 10.f));
        ImGui::PushStyleColor(ImGuiCol_MenuBarBg, menuBg);
        if (ImGui::BeginMainMenuBar()) {
            height = ImGui::GetWindowHeight();
            if (logoTexture) {
                const ImVec2 barPos = ImGui::GetWindowPos();
                const float logoSize = ImGui::GetTextLineHeight();
                ImGui::GetWindowDrawList()->AddImage(static_cast<ImTextureID>(logoTexture),
                    ImVec2(barPos.x + 10.f, barPos.y + (height - logoSize) * .5f),
                    ImVec2(barPos.x + 10.f + logoSize, barPos.y + (height + logoSize) * .5f));
            }
            ImGui::SetCursorPosX(logoTexture ? ImGui::GetTextLineHeight() + 18.f : 12.f);
            ImGui::TextColored(accent, "LIMA");
            ImGui::SameLine(0.f, 16.f);
            if (ImGui::BeginMenu("Representation")) {
                ColoringMenu(settings, commands);
                ImGui::Separator();
                if (ImGui::MenuItem("Show solvents", nullptr, &settings.showSolvents))
                    commands.push_back(Overlay::SolventVisibility{ settings.showSolvents });
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Camera")) {
				CameraMenu(commands, simulationId);
                if (!enableConsole) {
					ImGui::Separator();
                    if (ImGui::MenuItem("Single", nullptr, !tiled))
                        commands.push_back(Overlay::SetTiled{ false });
                    if (ImGui::MenuItem("Tiles", nullptr, tiled))
                        commands.push_back(Overlay::SetTiled{ true });
                }
                ImGui::EndMenu();
            }
            const auto fpsText = Lima::Format("{} fps", fps);
            const float fpsX = ImGui::GetWindowWidth() - ImGui::CalcTextSize(fpsText.c_str()).x - 14.f;
            if (fpsX > ImGui::GetCursorPosX() + 20.f) {
                ImGui::SetCursorPosX(fpsX);
                ImGui::TextColored(mutedText, "%s", fpsText.c_str());
            }
            ImGui::EndMainMenuBar();
            ImGui::GetForegroundDrawList()->AddLine(
                ImVec2(0.f, height - 1.f),
                ImVec2(ImGui::GetIO().DisplaySize.x, height - 1.f),
                ImGui::GetColorU32(panelBorder)
            );
        }
        ImGui::PopStyleColor();
        ImGui::PopStyleVar();
        return height;
    }

    float DrawSimulationTabs(const std::vector<SimulationTab>& tabs,
        std::deque<Overlay::Command>& commands, float top)
    {
        if (tabs.size() < 2)
            return top + 12.f;
        const float width = ImGui::GetIO().DisplaySize.x;
        const float tabHeight = ImGui::GetTextLineHeight() + 8.f;
        float totalWidth = 0.f;
        for (const auto& tab : tabs)
            totalWidth += ImGui::CalcTextSize(tab.label.c_str()).x + 24.f;
        const bool overflow = totalWidth > width;
        const float height = tabHeight + (overflow ? ImGui::GetStyle().ScrollbarSize : 0.f);

        ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0.f, 0.f));
        ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.f);
        ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.f);
        ImGui::PushStyleVar(ImGuiStyleVar_WindowMinSize, ImVec2(1.f, 1.f));
        ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(12.f, 4.f));
        ImGui::PushStyleVar(ImGuiStyleVar_FrameRounding, 0.f);
        ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(0.f, 0.f));
        ImGui::PushStyleColor(ImGuiCol_WindowBg, tabBg);
        BeginPanel("##SimulationTabs", ImVec2(0.f, top), ImVec2(width, height),
            overflow ? ImGuiWindowFlags_HorizontalScrollbar : ImGuiWindowFlags_NoScrollbar);
        for (size_t i = 0; i < tabs.size(); ++i) {
            if (i > 0)
                ImGui::SameLine();
            const auto& tab = tabs[i];
            ImGui::PushID(static_cast<int>(i));
            ImGui::PushStyleColor(ImGuiCol_Button, tab.active ? selectedBg : ImVec4(0.f, 0.f, 0.f, 0.f));
            ImGui::PushStyleColor(ImGuiCol_Text, tab.active ? textColor : mutedText);
            const float tabWidth = ImGui::CalcTextSize(tab.label.c_str()).x + 24.f;
            if (ImGui::Button(tab.label.c_str(), ImVec2(tabWidth, tabHeight)) && !tab.active)
                commands.push_back(Overlay::SelectSimulation{ tab.simulationId });
            Tooltip(tab.completed ? "Completed simulation" : "Simulation");
            const auto min = ImGui::GetItemRectMin();
            const auto max = ImGui::GetItemRectMax();
            if (tab.active)
                ImGui::GetWindowDrawList()->AddLine(ImVec2(min.x, min.y + 1.f),
                    ImVec2(max.x, min.y + 1.f), ImGui::GetColorU32(accent), 2.f);
            ImGui::GetWindowDrawList()->AddLine(ImVec2(max.x - 1.f, min.y + 6.f),
                ImVec2(max.x - 1.f, max.y - 6.f), ImGui::GetColorU32(panelBorder));
            ImGui::PopStyleColor(2);
            ImGui::PopID();
        }
        ImGui::End();
        ImGui::PopStyleColor();
        ImGui::PopStyleVar(7);
        return top + height + 12.f;
    }
}

Overlay::Overlay(GLFWwindow* window, const std::filesystem::path& limaDir)
{
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();

    ImGuiIO& io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
    io.IniFilename = nullptr;

    float contentScaleX = 1.0f;
    float contentScaleY = 1.0f;
    glfwGetWindowContentScale(window, &contentScaleX, &contentScaleY);
    const float contentScale = std::max(contentScaleX, contentScaleY);

    // Better default choice than Roboto for this kind of UI.
    // Put Inter-Medium.ttf in resources/ui if you have it.
    if (std::filesystem::exists(limaDir / "resources" / "ui" / "Inter-Medium.ttf")) {
        io.Fonts->AddFontFromFileTTF(
            (limaDir / "resources" / "ui" / "Inter-Medium.ttf").string().c_str(),
            22.0f * contentScale
        );
    }
    else {
        io.Fonts->AddFontFromFileTTF(
            (limaDir / "resources" / "ui" / "Roboto-Medium.ttf").string().c_str(),
            22.0f * contentScale
        );
    }

    PushOverlayTheme();
    ImGui::GetStyle().FontScaleMain = 1.0f / contentScale;

    int logoWidth = 0;
    int logoHeight = 0;
    int logoChannels = 0;
    unsigned char* logoPixels = stbi_load(
        (limaDir / "resources" / "logo" / "Lima_Symbol_64x64.png").string().c_str(),
        &logoWidth, &logoHeight, &logoChannels, 4);
    if (logoPixels) {
        glGenTextures(1, &logoTexture);
        glBindTexture(GL_TEXTURE_2D, logoTexture);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, logoWidth, logoHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, logoPixels);
        stbi_image_free(logoPixels);
    }

    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 430");
}

Overlay::~Overlay()
{
    if (logoTexture)
        glDeleteTextures(1, &logoTexture);
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
}


void Overlay::HandleConsole()
{
    const auto displaySize = ImGui::GetIO().DisplaySize;
    const float width = std::min(920.f, displaySize.x - margin * 2.f);
    BeginPanel("##CommandConsole",
        ImVec2((displaySize.x - width) * .5f, displaySize.y - margin - consoleHeight),
        ImVec2(width, consoleHeight));
    SectionLabel("COMMAND CONSOLE");
    ImGui::BeginChild("##ConsoleHistory", ImVec2(0.f, -ImGui::GetFrameHeightWithSpacing()));
    for (const auto& line : consoleLines)
        ImGui::TextWrapped("%s", line.c_str());
    if (scrollConsoleToBottom) {
        ImGui::SetScrollHereY(1.f);
        scrollConsoleToBottom = false;
    }
    ImGui::EndChild();
    ImGui::SetNextItemWidth(-1.f);
    if (ImGui::InputTextWithHint("##Command", "Enter a command...", consoleInput.data(), consoleInput.size(),
        ImGuiInputTextFlags_EnterReturnsTrue)) {
        if (consoleInput[0] != '\0') {
            const std::string command{ consoleInput.data() };
            consoleLines.push_back("> " + command);
            if (consoleLines.size() > 100)
                consoleLines.pop_front();
            submittedCommands.push_back(SubmittedCmd{ command });
            consoleInput[0] = '\0';
            scrollConsoleToBottom = true;
        }
        ImGui::SetKeyboardFocusHere(-1);
    }
    ImGui::End();
}

void Overlay::HandleContextMenu(RenderSettings& settings, std::optional<glm::dvec2> pos, SimulationId simulationId)
{
    if (pos) {
        ImGui::SetNextWindowPos(ImVec2(static_cast<float>(pos->x), static_cast<float>(pos->y)));
        ImGui::OpenPopup("##ViewportMenu");
    }
    if (ImGui::BeginPopup("##ViewportMenu")) {
        SectionLabel("APPEARANCE");
    ColoringMenu(settings, submittedCommands);
        ImGui::Separator();
        CameraMenu(submittedCommands, simulationId);
        ImGui::EndPopup();
    }
}

float Overlay::BeginFrame(RenderSettings& settings, int fps, const std::vector<SimulationTab>& tabs,
    bool tiled, SimulationId simulationId)
{
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
    const float menuHeight = DrawMenuBar(settings, submittedCommands, fps, logoTexture,
        tiled, enableConsole, simulationId);
    return tiled || tabs.empty() ? menuHeight : DrawSimulationTabs(tabs, submittedCommands, menuHeight);
}

void Overlay::DrawTile(SimulationId simulationId, const RenderContext& context,
    const RenderViewport& viewport, bool tiled)
{
    const float scale = tiled ? std::clamp(static_cast<float>(std::min(viewport.size.x, viewport.size.y)) / 600.f, .45f, .85f) : 1.f;
    const float top = static_cast<float>(viewport.origin.y) + 12.f * scale;
    const float bottom = static_cast<float>(viewport.origin.y + viewport.size.y) - 12.f * scale
        - (enableConsole ? consoleHeight + margin : 0.f);
    cardLayout = { top, top, bottom, static_cast<float>(viewport.origin.x),
        static_cast<float>(viewport.size.x), tiled, simulationId, scale };
    ImGui::PushFont(nullptr, ImGui::GetStyle().FontSizeBase * scale);
    if (tiled) {
        const auto label = context.label.empty() ? Lima::Format("Simulation {}", simulationId + 1) : context.label;
        const auto name = Lima::Format("##TileTitle{}", simulationId);
        const float titleHeight = ImGui::GetTextLineHeight() + 6.f * scale;
        ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(8.f * scale, 3.f * scale));
        ImGui::PushStyleVar(ImGuiStyleVar_WindowMinSize, ImVec2(1.f, 1.f));
        ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.f);
        BeginPanel(name.c_str(), ImVec2(static_cast<float>(viewport.origin.x), static_cast<float>(viewport.origin.y)),
            ImVec2(static_cast<float>(viewport.size.x), titleHeight), ImGuiWindowFlags_NoInputs);
        ImGui::TextUnformatted(label.c_str());
        if (context.completed) {
            ImGui::SameLine();
            ImGui::TextDisabled("(completed)");
        }
        ImGui::End();
        ImGui::PopStyleVar(3);
        cardLayout.leftY = cardLayout.rightY = top + titleHeight;
    }
    if (const auto* task = std::get_if<std::unique_ptr<Rendering::AtomRenderTask>>(&context.currentRenderTask))
        DrawTelemetry((*task)->simStatus);
    DrawMoleculeInfo(context.selectedMolecule);
    ImGui::PopFont();
}


void Overlay::EndFrame(RenderSettings& settings, std::optional<glm::dvec2> rightClickedPos,
    std::optional<SimulationId> popupSimulationId, bool spinnerVisible)
{
    if (enableConsole)
        HandleConsole();
    if (popupSimulationId)
        HandleContextMenu(settings, rightClickedPos, *popupSimulationId);
    if (spinnerVisible) {
        const auto displaySize = ImGui::GetIO().DisplaySize;
        const ImVec2 center(displaySize.x - 35.f, 70.f);
        auto* drawList = ImGui::GetForegroundDrawList();
        const float angle = static_cast<float>(ImGui::GetTime()) * 4.f;
        drawList->PathArcTo(center, 9.f, angle, angle + 4.7f, 24);
        drawList->PathStroke(ImGui::GetColorU32(accent), 0, 2.f);
    }
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
