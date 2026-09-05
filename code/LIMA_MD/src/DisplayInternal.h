#pragma once

#include "Display.h"
#include "RenderCommons.h"

#include <chrono>
#include <deque>
#include <filesystem>
#include <optional>
#include <string>
#include <variant>
#include <vector>
#include <array>

#include <glm.hpp>

class DrawTrianglesShader;
class GLFWwindow;

class FPS {
	std::array<std::chrono::high_resolution_clock::time_point, 32> prevTimepoints;
	int head = 0;
public:
	FPS();
	void NewFrame();
	int GetFps() const;
};

class Camera {
	Float3 center;
	Float3 boxSize;
	float dist = -2.f;
	float yaw = 0;
	float pitch = 0;
	float aspectRatio = 1.f;

public:
	Camera(Float3 boxSize);
	void Update(float deltaYaw, float deltaPitch, float deltaDist);
	void Update(Float3 boxSize);
	void UpdateViewport(glm::ivec2 viewportSize);
	void Reset();

	glm::mat4 View() const;
	glm::mat4 Projection() const;
	glm::mat4 ViewProjection() const;
};

struct RenderSettings {
	bool showSolvents = true;
	ColoringMethod coloringMethod{};
	bool hasForceData = false;
	bool hasBackbone = false;
};

class Overlay {
public:
	struct SubmittedCmd { std::string cmd{}; };
	struct SolventVisibility { bool visible = true; };
	struct ResetCamera {};
	struct RevolveCamera {};
	using Command = std::variant<SubmittedCmd, ColoringMethod, SolventVisibility, ResetCamera, RevolveCamera>;
private:
	bool didDrawThisFrame = false;

	void HandleConsole();
	void HandleContextMenu(RenderSettings& renderSettings, std::optional<glm::dvec2> rightClickedPos);

public:
	std::deque<Command> submittedCommands;
	bool enableConsole = false;

	Overlay(GLFWwindow*, const std::filesystem::path& limadir);
	~Overlay();

	void Draw(RenderSettings&, const SimStatus&, int fps,
		std::optional<glm::dvec2> rightClickedPos, bool spinnerVisible);
	void Render();
};

struct Arrow {
	glm::vec3 direction = glm::vec3(1.f, 0.f, 0.f);
	std::vector<Vertex> vertices;
	glm::vec4 color;
	int uniqueId;
	Arrow(glm::vec3 direction, glm::vec4 color, int uniqueId);
	void Draw(DrawTrianglesShader*, const glm::mat4& MVP, const glm::vec3& position, float scale = 1.f) const;
};

struct Ring {
	glm::vec3 normal = glm::vec3(0.f, 0.f, 1.f);
	std::vector<Vertex> vertices;
	glm::vec4 color;
	int uniqueId;
	Ring(glm::vec3 normal, glm::vec4 color, int uniqueId);
	void Draw(DrawTrianglesShader* shader, const glm::mat4& VP, const glm::vec3& position, float scale = 1.f) const;
};

struct TransformGizmo {
	glm::vec3 position{};
	int idOfAtomAttachedTo = -1;

	std::optional<int> activeAxis = std::nullopt;
	enum GizmoMode { Translate, Rotate } activeMode = GizmoMode::Translate;

	std::optional<glm::vec3> pullForce;
	std::optional<glm::vec3> rotateForce;

	glm::vec3 dragStartPosition{};
	glm::vec2 dragStartMousePos{};
	glm::vec3 dragStartRotateVector{};

	Arrow arrowX{ glm::vec3(1.f, 0.f, 0.f), glm::vec4(1.f, 0.f, 0.f, 1.f), static_cast<int>(UniqueRenderElementIds::gizmoArrowX) };
	Arrow arrowY{ glm::vec3(0.f, 1.f, 0.f), glm::vec4(0.f, 1.f, 0.f, 1.f), static_cast<int>(UniqueRenderElementIds::gizmoArrowY) };
	Arrow arrowZ{ glm::vec3(0.f, 0.f, 1.f), glm::vec4(0.f, 0.f, 1.f, 1.f), static_cast<int>(UniqueRenderElementIds::gizmoArrowZ) };

	Ring ringX{ glm::vec3(1.f, 0.f, 0.f), glm::vec4(1.f, 0.25f, 0.25f, 1.f), static_cast<int>(UniqueRenderElementIds::gizmoRotateX) };
	Ring ringY{ glm::vec3(0.f, 1.f, 0.f), glm::vec4(0.25f, 1.f, 0.25f, 1.f), static_cast<int>(UniqueRenderElementIds::gizmoRotateY) };
	Ring ringZ{ glm::vec3(0.f, 0.f, 1.f), glm::vec4(0.25f, 0.25f, 1.f, 1.f), static_cast<int>(UniqueRenderElementIds::gizmoRotateZ) };

	void Draw(DrawTrianglesShader* shader, const glm::mat4& VP) const;
	void SetActiveAxis(int selectedObjectId);
	void BeginDragging(glm::vec2 mousePos, const Camera& camera, glm::vec2 windowSize);
	void UpdateDraggingForce(glm::vec2 mousePos, const Camera& camera, glm::vec2 windowSize);
};
