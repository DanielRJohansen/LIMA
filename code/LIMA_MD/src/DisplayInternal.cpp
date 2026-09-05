#include "DisplayInternal.h"

#include "Utilities.h"

#include <algorithm>
#include <cmath>

#include <gtc/matrix_transform.hpp>

Camera::Camera(Float3 boxSize)
	: center(boxSize / 2.f)
	, boxSize(boxSize)
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
	this->boxSize = boxSize;
	if (center != boxSize / 2.f) {
		const float currentAspectRatio = aspectRatio;
		*this = Camera(boxSize);
		aspectRatio = currentAspectRatio;
	}
}

void Camera::Reset()
{
	const float currentAspectRatio = aspectRatio;
	*this = Camera(boxSize);
	aspectRatio = currentAspectRatio;
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




Arrow::Arrow(glm::vec3 direction, glm::vec4 color, int id) : direction(glm::normalize(direction)), color(color), uniqueId(id) {
	constexpr int radialSegments = 64;
	constexpr float totalLength = 1.0f;
	constexpr float shaftLength = 0.78f * totalLength;
	constexpr float shaftRadius = 0.035f * totalLength;
	constexpr float headLength = totalLength - shaftLength;
	constexpr float headRadius = 0.09f * totalLength;

	const glm::vec3 axis(0.f, 0.f, 1.f);
	const glm::vec3 u(1.f, 0.f, 0.f);
	const glm::vec3 v(0.f, 1.f, 0.f);

	const glm::vec3 shaftStart(0.f, 0.f, 0.f);
	const glm::vec3 shaftEnd(0.f, 0.f, shaftLength);
	const glm::vec3 coneBase(0.f, 0.f, shaftLength);
	const glm::vec3 apex(0.f, 0.f, totalLength);

	vertices.reserve(radialSegments * 12);

	auto AppendTriangle = [&](const glm::vec3& a, const glm::vec3& b, const glm::vec3& c) {
		const glm::vec3 normal = glm::normalize(glm::cross(b - a, c - a));
		vertices.push_back({ a, normal });
		vertices.push_back({ b, normal });
		vertices.push_back({ c, normal });
		};

	for (int i = 0; i < radialSegments; ++i) {
		const float a0 = 2.f * 3.1415f * static_cast<float>(i) / static_cast<float>(radialSegments);
		const float a1 = 2.f * 3.1415f * static_cast<float>(i + 1) / static_cast<float>(radialSegments);

		const glm::vec3 r0s = std::cos(a0) * u * shaftRadius + std::sin(a0) * v * shaftRadius;
		const glm::vec3 r1s = std::cos(a1) * u * shaftRadius + std::sin(a1) * v * shaftRadius;

		const glm::vec3 p0 = shaftStart + r0s;
		const glm::vec3 p1 = shaftStart + r1s;
		const glm::vec3 q0 = shaftEnd + r0s;
		const glm::vec3 q1 = shaftEnd + r1s;

		AppendTriangle(p0, p1, q1);
		AppendTriangle(p0, q1, q0);

		const glm::vec3 r0c = std::cos(a0) * u * headRadius + std::sin(a0) * v * headRadius;
		const glm::vec3 r1c = std::cos(a1) * u * headRadius + std::sin(a1) * v * headRadius;

		const glm::vec3 c0 = coneBase + r0c;
		const glm::vec3 c1 = coneBase + r1c;

		AppendTriangle(c0, c1, apex);
		AppendTriangle(coneBase, c1, c0);
	}

	AppendTriangle(shaftStart, shaftStart + glm::vec3(shaftRadius, 0.f, 0.f), shaftStart + glm::vec3(0.f, shaftRadius, 0.f));
}




Ring::Ring(glm::vec3 normal, glm::vec4 color, int id)
	: normal(glm::normalize(normal)), color(color), uniqueId(id)
{
	constexpr int majorSegments = 96;
	constexpr int minorSegments = 12;
	constexpr float majorRadius = 1.45f;
	constexpr float tubeRadius = 0.020f;

	vertices.reserve(majorSegments * minorSegments * 6);

	auto TorusPoint = [&](float u, float v) {
		const float cu = std::cos(u);
		const float su = std::sin(u);
		const float cv = std::cos(v);
		const float sv = std::sin(v);

		const float r = majorRadius + tubeRadius * cv;
		return glm::vec3(r * cu, r * su, tubeRadius * sv);
		};

	auto TorusNormal = [&](float u, float v) {
		const float cu = std::cos(u);
		const float su = std::sin(u);
		const float cv = std::cos(v);
		const float sv = std::sin(v);

		return glm::normalize(glm::vec3(cv * cu, cv * su, sv));
		};

	auto AppendTri = [&](const glm::vec3& a, const glm::vec3& na,
		const glm::vec3& b, const glm::vec3& nb,
		const glm::vec3& c, const glm::vec3& nc)
		{
			vertices.push_back({ a, na });
			vertices.push_back({ b, nb });
			vertices.push_back({ c, nc });
		};

	for (int i = 0; i < majorSegments; ++i) {
		const float u0 = 2.f * PI * static_cast<float>(i) / static_cast<float>(majorSegments);
		const float u1 = 2.f * PI * static_cast<float>(i + 1) / static_cast<float>(majorSegments);

		for (int j = 0; j < minorSegments; ++j) {
			const float v0 = 2.f * PI * static_cast<float>(j) / static_cast<float>(minorSegments);
			const float v1 = 2.f * PI * static_cast<float>(j + 1) / static_cast<float>(minorSegments);

			const glm::vec3 p00 = TorusPoint(u0, v0);
			const glm::vec3 p10 = TorusPoint(u1, v0);
			const glm::vec3 p11 = TorusPoint(u1, v1);
			const glm::vec3 p01 = TorusPoint(u0, v1);

			const glm::vec3 n00 = TorusNormal(u0, v0);
			const glm::vec3 n10 = TorusNormal(u1, v0);
			const glm::vec3 n11 = TorusNormal(u1, v1);
			const glm::vec3 n01 = TorusNormal(u0, v1);

			AppendTri(p00, n00, p10, n10, p11, n11);
			AppendTri(p00, n00, p11, n11, p01, n01);
		}
	}
}


void TransformGizmo::Draw(DrawTrianglesShader* shader, const glm::mat4& VP) const {
	const int axis = activeAxis.value_or(-1);

	const float translateScaleX = activeMode == GizmoMode::Translate && axis == 0 ? 2.2f : 2.f;
	const float translateScaleY = activeMode == GizmoMode::Translate && axis == 1 ? 2.2f : 2.f;
	const float translateScaleZ = activeMode == GizmoMode::Translate && axis == 2 ? 2.2f : 2.f;

	const float rotateScaleX = activeMode == GizmoMode::Rotate && axis == 0 ? 2.2f : 2.f;
	const float rotateScaleY = activeMode == GizmoMode::Rotate && axis == 1 ? 2.2f : 2.f;
	const float rotateScaleZ = activeMode == GizmoMode::Rotate && axis == 2 ? 2.2f : 2.f;

	arrowX.Draw(shader, VP, position, translateScaleX);
	arrowY.Draw(shader, VP, position, translateScaleY);
	arrowZ.Draw(shader, VP, position, translateScaleZ);

	ringX.Draw(shader, VP, position, rotateScaleX);
	ringY.Draw(shader, VP, position, rotateScaleY);
	ringZ.Draw(shader, VP, position, rotateScaleZ);
}