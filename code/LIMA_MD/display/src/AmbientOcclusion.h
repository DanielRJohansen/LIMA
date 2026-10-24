#pragma once

#include "GLShader.h"

#include <array>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <string>

// Screen-space ambient occlusion: darkens crevices and contacts, which is what gives packed molecules their depth.
// The scene is drawn into an offscreen color and depth texture. A first pass reconstructs view-space positions and
// normals from depth and counts how much of a hemisphere around each pixel is inside geometry. A second pass blurs
// that, without bleeding across depth edges, and writes the darkened color to the window.
// Works from depth alone, so it applies to everything drawn: atoms, lipid tubes and the backbone cartoon.
class AmbientOcclusion {
	static constexpr int kernelSize = 16;
	static constexpr float radius = 0.5f;	// [nm] Reach of the occlusion, about the size of an atom
	// Fraction of the occlusion applied. LIMA_AO_STRENGTH overrides it, for tuning and before/after comparisons
	static float Strength() {
		static const float strength = std::getenv("LIMA_AO_STRENGTH") ? std::stof(std::getenv("LIMA_AO_STRENGTH")) : 0.9f;
		return strength;
	}

	static constexpr const char* fullscreenVertexSource = R"(
#version 430 core
void main() {
    // One triangle that covers the viewport
    const vec2 corner = vec2((gl_VertexID << 1) & 2, gl_VertexID & 2);
    gl_Position = vec4(corner * 2.0 - 1.0, 0.0, 1.0);
}
)";

	static constexpr const char* occlusionFragmentSource = R"(
#version 430 core

uniform sampler2D Depth;
uniform mat4 Proj;
uniform mat4 InvProj;
uniform ivec4 Viewport;	// {x, y, width, height} in framebuffer pixels
uniform vec3 Kernel[16];
uniform float Radius;

layout(location = 0) out float Occlusion;

vec3 ViewPosition(ivec2 pixel) {
    const float depth = texelFetch(Depth, pixel, 0).r;
    const vec2 ndc = (vec2(pixel - Viewport.xy) + 0.5) / vec2(Viewport.zw) * 2.0 - 1.0;
    const vec4 view = InvProj * vec4(ndc, depth * 2.0 - 1.0, 1.0);
    return view.xyz / view.w;
}

// Per-pixel pseudo random angle, so the kernel's banding turns into noise that the blur removes
float Hash(vec2 p) {
    return fract(sin(dot(p, vec2(12.9898, 78.233))) * 43758.5453);
}

void main() {
    const ivec2 pixel = ivec2(gl_FragCoord.xy);
    if (texelFetch(Depth, pixel, 0).r >= 1.0) {
        Occlusion = 0.0;	// Background
        return;
    }
    const vec3 position = ViewPosition(pixel);
    // The normal from the depth of the neighbors, taking the smaller step on each axis to stay off silhouettes
    const vec3 right = ViewPosition(pixel + ivec2(1, 0)) - position;
    const vec3 left = position - ViewPosition(pixel - ivec2(1, 0));
    const vec3 up = ViewPosition(pixel + ivec2(0, 1)) - position;
    const vec3 down = position - ViewPosition(pixel - ivec2(0, 1));
    const vec3 dx = abs(right.z) < abs(left.z) ? right : left;
    const vec3 dy = abs(up.z) < abs(down.z) ? up : down;
    const vec3 normal = normalize(cross(dx, dy));

    const float angle = Hash(gl_FragCoord.xy) * 6.2831853;
    const vec3 random = vec3(cos(angle), sin(angle), 0.0);
    const vec3 tangent = normalize(random - normal * dot(random, normal));
    const mat3 tbn = mat3(tangent, cross(normal, tangent), normal);

    float occlusion = 0.0;
    for (int i = 0; i < 16; ++i) {
        const vec3 samplePosition = position + tbn * Kernel[i] * Radius;
        const vec4 clip = Proj * vec4(samplePosition, 1.0);
        const vec2 ndc = clip.xy / clip.w;
        const ivec2 samplePixel = Viewport.xy + ivec2((ndc * 0.5 + 0.5) * vec2(Viewport.zw));
        if (any(lessThan(samplePixel, Viewport.xy)) || any(greaterThanEqual(samplePixel, Viewport.xy + Viewport.zw)))
            continue;
        const float sceneZ = ViewPosition(samplePixel).z;
        // Occluded when the scene is in front of the sample, weighted down when it is far in front (another molecule)
        const float rangeWeight = smoothstep(0.0, 1.0, Radius / abs(position.z - sceneZ));
        occlusion += (sceneZ >= samplePosition.z + 0.02 ? 1.0 : 0.0) * rangeWeight;
    }
    Occlusion = occlusion / 16.0;
}
)";

	static constexpr const char* compositeFragmentSource = R"(
#version 430 core

uniform sampler2D Color;
uniform sampler2D Depth;
uniform sampler2D OcclusionMap;
uniform ivec4 Viewport;
uniform float Strength;

layout(location = 0) out vec4 FragColor;

void main() {
    const ivec2 pixel = ivec2(gl_FragCoord.xy);
    const vec4 color = texelFetch(Color, pixel, 0);
    const float depth = texelFetch(Depth, pixel, 0).r;
    if (depth >= 1.0) {
        FragColor = vec4(color.rgb, 1.0);
        return;
    }

    // 4x4 blur that skips samples across a depth edge, so occlusion stays inside its molecule
    float sum = 0.0;
    float weight = 0.0;
    for (int y = -2; y < 2; ++y) {
        for (int x = -2; x < 2; ++x) {
            const ivec2 p = clamp(pixel + ivec2(x, y), Viewport.xy, Viewport.xy + Viewport.zw - 1);
            const float w = abs(texelFetch(Depth, p, 0).r - depth) < 0.002 ? 1.0 : 0.0;
            sum += texelFetch(OcclusionMap, p, 0).r * w;
            weight += w;
        }
    }
    const float occlusion = sum / max(weight, 1.0);
    FragColor = vec4(color.rgb * (1.0 - Strength * occlusion), 1.0);
}
)";

	Shader occlusionShader{ fullscreenVertexSource, occlusionFragmentSource };
	Shader compositeShader{ fullscreenVertexSource, compositeFragmentSource };
	std::array<glm::vec3, kernelSize> kernel{};

	GLuint vao = 0;	// Core profile needs one bound to draw, even without vertex attributes
	GLuint sceneFramebuffer = 0;
	GLuint sceneColor = 0;
	GLuint sceneDepth = 0;
	GLuint occlusionFramebuffer = 0;
	GLuint occlusionTexture = 0;
	glm::ivec2 size{};

	static GLuint MakeTexture(GLenum internalFormat, GLenum format, GLenum type, glm::ivec2 size) {
		GLuint texture = 0;
		glGenTextures(1, &texture);
		glBindTexture(GL_TEXTURE_2D, texture);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexImage2D(GL_TEXTURE_2D, 0, internalFormat, size.x, size.y, 0, format, type, nullptr);
		return texture;
	}

	void Destroy() {
		const GLuint textures[] = { sceneColor, sceneDepth, occlusionTexture };
		glDeleteTextures(3, textures);
		const GLuint framebuffers[] = { sceneFramebuffer, occlusionFramebuffer };
		glDeleteFramebuffers(2, framebuffers);
		sceneColor = sceneDepth = occlusionTexture = sceneFramebuffer = occlusionFramebuffer = 0;
	}

	static void BindTexture(Shader& shader, const char* name, int unit, GLuint texture) {
		glActiveTexture(GL_TEXTURE0 + unit);
		glBindTexture(GL_TEXTURE_2D, texture);
		shader.SetUniformI(name, unit);
	}

	static void SetViewport(Shader& shader, glm::ivec4 viewport) {
		glUniform4i(glGetUniformLocation(shader.programID, "Viewport"), viewport.x, viewport.y, viewport.z, viewport.w);
	}

public:
	AmbientOcclusion() {
		glGenVertexArrays(1, &vao);
		// Samples in a hemisphere around +z, denser near the center, as close occluders matter most
		for (int i = 0; i < kernelSize; ++i) {
			const float golden = 2.39996323f * i;
			const float z = 1.f - (i + 0.5f) / kernelSize;
			const float r = std::sqrt(1.f - z * z);
			const float scale = 0.1f + 0.9f * std::pow((i + 1.f) / kernelSize, 2.f);
			kernel[i] = glm::vec3(r * std::cos(golden), r * std::sin(golden), z) * scale;
		}
	}
	~AmbientOcclusion() {
		Destroy();
		glDeleteVertexArrays(1, &vao);
	}

	void Resize(glm::ivec2 newSize) {
		if (newSize.x <= 0 || newSize.y <= 0 || (sceneFramebuffer && newSize == size))
			return;
		size = newSize;
		Destroy();

		sceneColor = MakeTexture(GL_RGBA8, GL_RGBA, GL_UNSIGNED_BYTE, size);
		sceneDepth = MakeTexture(GL_DEPTH_COMPONENT24, GL_DEPTH_COMPONENT, GL_UNSIGNED_INT, size);
		glGenFramebuffers(1, &sceneFramebuffer);
		glBindFramebuffer(GL_FRAMEBUFFER, sceneFramebuffer);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, sceneColor, 0);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, sceneDepth, 0);
		if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
			throw std::runtime_error("Ambient occlusion scene framebuffer incomplete");

		occlusionTexture = MakeTexture(GL_R8, GL_RED, GL_UNSIGNED_BYTE, size);
		glGenFramebuffers(1, &occlusionFramebuffer);
		glBindFramebuffer(GL_FRAMEBUFFER, occlusionFramebuffer);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, occlusionTexture, 0);
		if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
			throw std::runtime_error("Ambient occlusion framebuffer incomplete");
		glBindFramebuffer(GL_FRAMEBUFFER, 0);
	}

	// Draw the scene after this, with the viewport and scissor of its tile already set
	void BindScene() {
		glBindFramebuffer(GL_FRAMEBUFFER, sceneFramebuffer);
	}

	// Writes the scene, darkened by its occlusion, to the window's framebuffer in the same viewport
	void Resolve(const glm::mat4& projection, glm::ivec4 viewport) {
		glDisable(GL_DEPTH_TEST);
		glDisable(GL_BLEND);	// Copies the already blended scene
		glBindVertexArray(vao);

		glBindFramebuffer(GL_FRAMEBUFFER, occlusionFramebuffer);
		occlusionShader.use();
		BindTexture(occlusionShader, "Depth", 0, sceneDepth);
		occlusionShader.SetUniformMat4("Proj", projection);
		occlusionShader.SetUniformMat4("InvProj", glm::inverse(projection));
		occlusionShader.SetUniform("Radius", radius);
		glUniform3fv(glGetUniformLocation(occlusionShader.programID, "Kernel"), kernelSize, &kernel[0].x);
		SetViewport(occlusionShader, viewport);
		glDrawArrays(GL_TRIANGLES, 0, 3);

		glBindFramebuffer(GL_FRAMEBUFFER, 0);
		compositeShader.use();
		BindTexture(compositeShader, "Color", 0, sceneColor);
		BindTexture(compositeShader, "Depth", 1, sceneDepth);
		BindTexture(compositeShader, "OcclusionMap", 2, occlusionTexture);
		compositeShader.SetUniform("Strength", Strength());
		SetViewport(compositeShader, viewport);
		glDrawArrays(GL_TRIANGLES, 0, 3);

		glActiveTexture(GL_TEXTURE0);
		glBindVertexArray(0);
		glUseProgram(0);
		glEnable(GL_DEPTH_TEST);
		glEnable(GL_BLEND);
	}
};
