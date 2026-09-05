#include "DisplayInternal.h"

#include "Utilities.h"

#include <algorithm>
#include <cmath>

#include <gtc/matrix_transform.hpp>

Camera::Camera(Float3 boxSize)
	: center(boxSize / 2.f)
	, dist(-2.f * boxSize.y)
{}

void Camera::Update(float deltaYaw, float deltaPitch, float deltaDist)
{
	yaw += deltaYaw;
	pitch += deltaPitch;
	dist += deltaDist + deltaDist * -std::min(dist, 0.f) * 0.5f;
}

void Camera::Update(Float3 boxSize)
{
	if (center != boxSize / 2.f) {
		const float currentAspectRatio = aspectRatio;
		*this = Camera(boxSize);
		aspectRatio = currentAspectRatio;
	}
}

void Camera::UpdateViewport(glm::ivec2 viewportSize)
{
	if (viewportSize.x > 0 && viewportSize.y > 0)
		aspectRatio = static_cast<float>(viewportSize.x) / static_cast<float>(viewportSize.y);
}

glm::mat4 Camera::View() const
{
	glm::mat4 view(1.f);
	view = glm::translate(view, glm::vec3(0.f, 0.f, dist));
	view = glm::rotate(view, -PI / 2.f, glm::vec3(1.f, 0.f, 0.f));
	view = glm::rotate(view, pitch, glm::vec3(1.f, 0.f, 0.f));
	view = glm::rotate(view, yaw, glm::vec3(0.f, 0.f, 1.f));
	view = glm::translate(view, ToVec3(-center));
	return view;
}

glm::mat4 Camera::Projection() const
{
	constexpr double fovY = 45.;
	constexpr double nearPlane = 0.1;
	constexpr double farPlane = 1000.;
	const double halfHeight = std::tan(glm::radians(fovY / 2.)) * nearPlane;
	const double halfWidth = halfHeight * aspectRatio;
	return glm::frustum(-halfWidth, halfWidth, -halfHeight, halfHeight, nearPlane, farPlane);
}

glm::mat4 Camera::ViewProjection() const
{
	return Projection() * View();
}

FPS::FPS()
{
	const auto now = std::chrono::high_resolution_clock::now();
	std::ranges::fill(prevTimepoints, now);
}

void FPS::NewFrame()
{
	head = (head + 1) % prevTimepoints.size();
	prevTimepoints[head] = std::chrono::high_resolution_clock::now();
}

int FPS::GetFps() const
{
	const int back = (head + 1) % prevTimepoints.size();
	const auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(prevTimepoints[head] - prevTimepoints[back]);
	const auto averageFrameTime = elapsed / (prevTimepoints.size() - 1);
	return static_cast<int>(1e9 / averageFrameTime.count());
}
