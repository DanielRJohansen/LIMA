// No pragma, once include in Display.cpp
#include <GL/glew.h>

#include "MoleculeHull.cuh"
#include "GLShader.h"
#include <glm.hpp>
#include <gtc/matrix_transform.hpp>
#include <gtc/type_ptr.hpp>
#include "RenderCommons.h"
#include "SSBO.h"

#include <cuda_gl_interop.h>

#include <mutex>
#include <condition_variable>
#include <optional>
#include <atomic>
#include <future>


class RenderTargetControl {
public:
    GLuint framebuffer = 0;
    GLuint colorTexture = 0;
    GLuint idTexture = 0;
    GLuint depthBuffer = 0;
    glm::ivec2 size{};

    RenderTargetControl() = default;
    ~RenderTargetControl() {
        Destroy();
    }

    struct ScopedDrawBinding {
        GLint prevDrawFbo = 0;
        GLint prevReadFbo = 0;
        GLint prevViewport[4]{};

        ScopedDrawBinding() {
            glGetIntegerv(GL_DRAW_FRAMEBUFFER_BINDING, &prevDrawFbo);
            glGetIntegerv(GL_READ_FRAMEBUFFER_BINDING, &prevReadFbo);
            glGetIntegerv(GL_VIEWPORT, prevViewport);
        }

        ~ScopedDrawBinding() {
            glBindFramebuffer(GL_DRAW_FRAMEBUFFER, prevDrawFbo);
            glBindFramebuffer(GL_READ_FRAMEBUFFER, prevReadFbo);
            glViewport(prevViewport[0], prevViewport[1], prevViewport[2], prevViewport[3]);
        }
    };

    void Destroy() {
        if (depthBuffer) { glDeleteRenderbuffers(1, &depthBuffer); depthBuffer = 0; }
        if (idTexture) { glDeleteTextures(1, &idTexture); idTexture = 0; }
        if (colorTexture) { glDeleteTextures(1, &colorTexture); colorTexture = 0; }
        if (framebuffer) { glDeleteFramebuffers(1, &framebuffer); framebuffer = 0; }
    }

    void Resize(glm::ivec2 newSize) {
		if (newSize.x <= 0 || newSize.y <= 0)
			return;
		if (framebuffer && size == newSize)
			return;

        size = newSize;
        Destroy();

        glGenFramebuffers(1, &framebuffer);
        glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);

        glGenTextures(1, &colorTexture);
        glBindTexture(GL_TEXTURE_2D, colorTexture);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, size.x, size.y, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, colorTexture, 0);

        glGenTextures(1, &idTexture);
        glBindTexture(GL_TEXTURE_2D, idTexture);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_R32I, size.x, size.y, 0, GL_RED_INTEGER, GL_INT, nullptr);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, idTexture, 0);

        glGenRenderbuffers(1, &depthBuffer);
        glBindRenderbuffer(GL_RENDERBUFFER, depthBuffer);
        glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, size.x, size.y);
        glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, depthBuffer);

        GLenum drawBuffers[] = { GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1 };
        glDrawBuffers(2, drawBuffers);

        if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
            glBindFramebuffer(GL_FRAMEBUFFER, 0);
            throw std::runtime_error("Picking framebuffer incomplete.");
        }

        glBindFramebuffer(GL_FRAMEBUFFER, 0);
    }

    [[nodiscard]]
    ScopedDrawBinding BindForDraw() {
        if (!framebuffer) {
            throw std::runtime_error("RenderTargetControl::BindForDraw: framebuffer not initialized.");
        }

        ScopedDrawBinding state{};

        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, framebuffer);
        glViewport(0, 0, size.x, size.y);

        return state;
    }

    void ClearForPicking() {
        if (!framebuffer) {
            throw std::runtime_error("RenderTargetControl::ClearForPicking: framebuffer not initialized.");
        }

        const GLfloat clearColor[4] = { 0.f, 0.f, 0.f, 0.f };
        glClearBufferfv(GL_COLOR, 0, clearColor);

        const GLint clearId[1] = { -1 };
        glClearBufferiv(GL_COLOR, 1, clearId);

        glClear(GL_DEPTH_BUFFER_BIT);
    }

    int ReadIdAtPixel(glm::ivec2 pixel) const {
        GLint prevReadFbo = 0;
        glGetIntegerv(GL_READ_FRAMEBUFFER_BINDING, &prevReadFbo);

        glBindFramebuffer(GL_READ_FRAMEBUFFER, framebuffer);
        glReadBuffer(GL_COLOR_ATTACHMENT1);

        int pixelValue = -1;
        glReadPixels(pixel.x,
            size.y - 1 - pixel.y,
            1, 1,
            GL_RED_INTEGER,
            GL_INT,
            &pixelValue);

        glBindFramebuffer(GL_READ_FRAMEBUFFER, prevReadFbo);
        return pixelValue;
    }
};

class DrawBackgroundGradientShader : public Shader {
    static constexpr const char* vertexShaderSource = R"(
        #version 430 core

        out vec2 uv;

        void main() {
            const vec2 positions[3] = vec2[3](
                vec2(-1.0, -1.0),
                vec2( 3.0, -1.0),
                vec2(-1.0,  3.0)
            );

            vec2 pos = positions[gl_VertexID];
            gl_Position = vec4(pos, 0.0, 1.0);
            uv = pos * 0.5 + 0.5;
        }
    )";

    static constexpr const char* fragmentShaderSource = R"(
        #version 430 core

        in vec2 uv;

        uniform vec4 colorBot;
        uniform vec4 colorTop;

        out vec4 FragColor;

        void main() {
            float t = clamp(uv.y, 0.0, 1.0);
            vec4 color = mix(colorBot, colorTop, t);
            FragColor = color;
        }
    )";

    GLuint VAO{};

public:
    DrawBackgroundGradientShader() : Shader(vertexShaderSource, fragmentShaderSource) {
        glGenVertexArrays(1, &VAO);
    }

    ~DrawBackgroundGradientShader() {
        glDeleteVertexArrays(1, &VAO);
    }

    void Draw(glm::vec4 colorBot, glm::vec4 colorTop) {
        use();

        const GLboolean depthTestWasEnabled = glIsEnabled(GL_DEPTH_TEST);
        GLint previousDepthMask = GL_TRUE;
        glGetIntegerv(GL_DEPTH_WRITEMASK, &previousDepthMask);

        glDisable(GL_DEPTH_TEST);
        glDepthMask(GL_FALSE);

        glBindVertexArray(VAO);
        SetUniform("colorBot", colorBot);
        SetUniform("colorTop", colorTop);
        glDrawArrays(GL_TRIANGLES, 0, 3);
        glBindVertexArray(0);

        glDepthMask(previousDepthMask);
        if (depthTestWasEnabled)
            glEnable(GL_DEPTH_TEST);

        glUseProgram(0);
    }
};

class DrawBoxOutlineShader : public Shader {
    static constexpr const char* vertexShaderSource = R"(
        #version 430 core
        layout(location = 0) in vec3 aPos;
        layout(location = 1) in vec3 aColor;

        out vec3 vertexColor;
        uniform mat4 MVP;
        uniform vec3 boxSize;

        void main() {
            gl_Position = MVP * vec4(aPos*boxSize, 1.0);
            vertexColor = aColor;
        }
        )";
    static constexpr const char* fragmentShaderSource = R"(
        #version 430 core
        in vec3 vertexColor;
        out vec4 FragColor;
        void main() {
            FragColor = vec4(vertexColor, 1.0);
        }
    )";

    GLuint VAO, VBO, EBO;

public:
    DrawBoxOutlineShader() : Shader(vertexShaderSource, fragmentShaderSource) {
        const float boxVertices[] = {
            // positions          // colors
            0.f, 0.f, 0.f,  0.2f, 0.2f, 0.8f,
            1.f, 0.f, 0.f,  0.2f, 0.2f, 0.8f,
            1.f, 1.f, 0.f,  0.2f, 0.2f, 0.8f,
            0.f, 1.f, 0.f,  0.2f, 0.2f, 0.8f,
            0.f, 0.f, 1.f,  0.4f, 0.4f, 0.4f,
            1.f, 0.f, 1.f,  0.4f, 0.4f, 0.4f,
            1.f, 1.f, 1.f,  0.8f, 0.4f, 0.4f,
            0.f, 1.f, 1.f,  0.8f, 0.4f, 0.4f,
        };

        const unsigned int boxIndices[] = {
            0, 1, 1, 2, 2, 3, 3, 0, // bottom
            4, 5, 5, 6, 6, 7, 7, 4, // top
            0, 4, 1, 5, 2, 6, 3, 7  // sides
        };

        glGenVertexArrays(1, &VAO);
        glGenBuffers(1, &VBO);
        glGenBuffers(1, &EBO);

        glBindVertexArray(VAO);

        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(boxVertices), boxVertices, GL_STATIC_DRAW);

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(boxIndices), boxIndices, GL_STATIC_DRAW);

        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
        glEnableVertexAttribArray(1);

        glBindVertexArray(0);
    }

    ~DrawBoxOutlineShader() {
		glDeleteVertexArrays(1, &VAO);
		glDeleteBuffers(1, &VBO);
		glDeleteBuffers(1, &EBO);
	}

    void Draw(const glm::mat4 MVP, Float3 boxSize) {
        use();

        glBindVertexArray(VAO);
        SetUniformMat4("MVP", MVP);
        SetUniformFloat3("boxSize", boxSize);
        glDrawElements(GL_LINES, 24, GL_UNSIGNED_INT, 0);
        glBindVertexArray(0);

        glUseProgram(0);
    }
};

class DrawFacetsShader : public Shader {
    static constexpr const char* hullVertexShaderSource = R"(
    #version 430 core

    struct Facet {
        float data[16]; // DONT CHANGE, it's OpenGL being annoying with offsets
        // vec3[3] vertices;
        // vec3 normal
        // float distFromOrigo
        // char[12] paddingBytes
    };

    layout(std430, binding = 0) buffer TriBuffer {
        Facet facets[];
    };

    uniform mat4 MVP;
    out vec4 fragColor;

    const vec3 lightDir = vec3(0.0, 0.0, -1.0); // Light coming from directly above
    const float colorAlpha = 0.5f;

    uniform vec3 boxSize;    
    uniform int drawMode; // 0 = FACES, 1 = EDGES
    uniform vec4 predecidedColor = vec4(-1.,-1.,-1.,-1.);  

    vec3 GenerateRandomColor(uint triIndex) {
		return vec3(
			fract(float(triIndex+1) * 0.6180339887498949), // Golden ratio conjugate
			fract(float(triIndex+1) * 0.7548776662466927), // Another irrational number
			fract(float(triIndex+1) * 0.514229)            // Fibonacci number
		);
	}

    void main() {
        uint triIndex = 0;
        uint vertexIndex = 0;
        
        if (drawMode == 0) {
            // For FACES mode, process vertices normally
            triIndex = gl_VertexID / 3;
            vertexIndex = gl_VertexID % 3;
        } else {
            // For EDGES mode, handle edge-specific vertex calculation
            triIndex = gl_VertexID / 6;
            uint edgeIndex = (gl_VertexID / 2) % 3;
            vertexIndex = (edgeIndex + gl_VertexID % 2) % 3;
        }
        

        // First set position     
        vec3 position = vec3(facets[triIndex].data[vertexIndex * 3], 
                             facets[triIndex].data[vertexIndex * 3 + 1], 
                             facets[triIndex].data[vertexIndex * 3 + 2]);// / boxSize - vec3(0.5f,0.5f,0.5f);
        gl_Position = MVP * vec4(position, 1.0);



        // Now set color
        vec3 triNormal = vec3(
                        facets[triIndex].data[3*3+0], 
                        facets[triIndex].data[3*3+1], 
                        facets[triIndex].data[3*3+2]);

        // Simulate area light by blending the normal with the light direction
        float brightness = clamp(
            dot(-triNormal, lightDir) * 0.5f + 0.5f, 
            0.1f,
            1.f);
        
        if (drawMode == 1) { // For EDGES mode, color the edges red  			
			fragColor = vec4(.8, 0.8, 0.8, colorAlpha);
		}
        else { // Generate a pseudo-random color based on the triangle index        
            fragColor = vec4(brightness * GenerateRandomColor(triIndex), 1.f);
        }
        
        if (predecidedColor.x != -1) {
			fragColor = predecidedColor;
		}
    }
)";


    static constexpr const char* hullFragmentShaderSource = R"(
    #version 430 core

    in vec4 fragColor;
    out vec4 color;

    void main() {
        color = fragColor;
    }
)";

    GLuint VAO;
    SSBO facetsBuffer{};

public:


    DrawFacetsShader() : Shader(hullVertexShaderSource, hullFragmentShaderSource) {
        // Generate and bind VAO
        glGenVertexArrays(1, &VAO);
        glBindVertexArray(VAO);
        glBindVertexArray(0);  // Unbind VAO
    }

    ~DrawFacetsShader() {
        glDeleteVertexArrays(1, &VAO);
    }

    void Draw(const glm::mat4& MVP, const Facet* facets_cudaMem, int numFacets, 
        FacetDrawMode mode, Float3 boxSize, std::optional<float4> color=std::nullopt)
    {
        use(); // Use the shader program

        // Update SSBO with CUDA device memory
        facetsBuffer.Bind(0);
        facetsBuffer.SetData_FromCuda(facets_cudaMem, numFacets * sizeof(Facet));

        // Unbind the buffer
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

        SetUniformFloat3("boxSize", boxSize);
        SetUniformMat4("MVP", MVP);
        SetUniformI("drawMode", mode == FACES ? 0 : 1);
        if (color.has_value()) {
			SetUniformFloat4("predecidedColor", color.value());
		}

        // Draw the triangles based on the selected mode
        glBindVertexArray(VAO);
        if (mode == FACES) {
            glDrawArrays(GL_TRIANGLES, 0, numFacets * 3);
        }
        else if (mode == EDGES) {
            glDrawArrays(GL_LINES, 0, numFacets * 6); // Each triangle has 3 edges, each edge has 2 vertices
        }
        glBindVertexArray(0);
        glUseProgram(0);  // Unbind shader program
    }
};

class DrawNormalsShader : public Shader {
    static constexpr const char* hullVertexShaderSource = R"(
    #version 430 core

    struct Facet {
        float data[16]; // DONT CHANGE, it's OpenGL being annoying with offsets
        // vec3[3] vertices;
        // vec3 normal
        // float distFromOrigo
        // char[12] paddingBytes
    };

    layout(std430, binding = 0) buffer TriBuffer {
        Facet facets[];
    };

    uniform mat4 MVP;
    out vec4 fragColor;

    uniform vec3 boxSize;

    vec3 GetFacetCenter(uint triIndex) {
		return (
            vec3(facets[triIndex].data[0], facets[triIndex].data[1], facets[triIndex].data[2])
            + vec3(facets[triIndex].data[3], facets[triIndex].data[4], facets[triIndex].data[5])
            + vec3(facets[triIndex].data[6], facets[triIndex].data[7], facets[triIndex].data[8])
		) / vec3(3.0f);
	}

    void main() {
        uint triIndex = gl_VertexID / 2; // Each normal is drawn with two vertices
        bool isTip = (gl_VertexID % 2 == 1); // Even index is base, odd index is tip
        


    
        vec3 position = GetFacetCenter(triIndex) / boxSize - vec3(0.5f, 0.5f, 0.5f);

        vec3 triNormal = vec3(facets[triIndex].data[3*3+0], 
                              facets[triIndex].data[3*3+1], 
                              facets[triIndex].data[3*3+2]);

        if (isTip) {
            position += triNormal * 0.02; // Move the tip 10% of the normal length away from the base
        }

        gl_Position = MVP * vec4(position, 1.0);

        // Set color based on the normal, could be replaced with a more advanced coloring
        fragColor = vec4(1,0,0,1);
    }
)";


    static constexpr const char* hullFragmentShaderSource = R"(
    #version 430 core

    in vec4 fragColor;
    out vec4 color;

    void main() {
        color = fragColor;
    }
)";

    GLuint VAO;
    SSBO facetsBuffer{};

public:
    DrawNormalsShader() : Shader(hullVertexShaderSource, hullFragmentShaderSource) {
        glGenVertexArrays(1, &VAO);
        glBindVertexArray(VAO);
        glBindVertexArray(0);
    }

    ~DrawNormalsShader() {
        glDeleteVertexArrays(1, &VAO);
    }

    void Draw(const glm::mat4& MVP, const Facet* facets_cudaMem, int numFacets, Float3 boxSize) {
        use();

        facetsBuffer.Bind(0);
        facetsBuffer.SetData_FromCuda(facets_cudaMem, numFacets * sizeof(Facet));

        //void* ptr = glMapBuffer(GL_SHADER_STORAGE_BUFFER, GL_WRITE_ONLY);
        //cudaMemcpy(ptr, facets_cudaMem, numFacets * sizeof(Facet), cudaMemcpyDeviceToHost);
        //glUnmapBuffer(GL_SHADER_STORAGE_BUFFER);

        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

        SetUniformFloat3("boxSize", boxSize);
        SetUniformMat4("MVP", MVP);

        glBindVertexArray(VAO);
        glDrawArrays(GL_LINES, 0, numFacets * 2); // Each normal is drawn with 2 vertices
        glBindVertexArray(0);
        glUseProgram(0);
    }
};

class DrawTrianglesShader final : public Shader {
private:
    static constexpr const char* vertexSource = R"(
        #version 430 core

        layout(location = 0) in vec3 inPosition;
        layout(location = 1) in vec3 inNormal;

        uniform mat4 MVP;
        uniform mat4 Model;

        out vec3 fragNormal;

        void main()
        {
            gl_Position = MVP * vec4(inPosition, 1.0);
            fragNormal = mat3(Model) * inNormal;
        }
    )";

    static constexpr const char* fragmentSource = R"(
        #version 430 core

        in vec3 fragNormal;

        uniform vec3 LightDir;
        uniform vec4 Color;
        uniform int ObjectId;

        layout(location = 0) out vec4 outColor;
        layout(location = 1) out int outObjectId;

        void main()
        {
            vec3 N = normalize(fragNormal);
            vec3 L = normalize(-LightDir);

            float diffuse = max(dot(N, L), 0.0);
            float ambient = 0.25;
            float lighting = ambient + diffuse * 0.75;

            outColor = vec4(Color.rgb * lighting, Color.a);
            outObjectId = ObjectId;
        }
    )";

    GLuint vao = 0;
    GLuint vbo = 0;

public:
    DrawTrianglesShader() : Shader(vertexSource, fragmentSource)
    {
        glCreateVertexArrays(1, &vao);
        glCreateBuffers(1, &vbo);

        glVertexArrayVertexBuffer(vao, 0, vbo, 0, sizeof(Vertex));

        glEnableVertexArrayAttrib(vao, 0);
        glEnableVertexArrayAttrib(vao, 1);

        glVertexArrayAttribFormat(vao, 0, 3, GL_FLOAT, GL_FALSE, offsetof(Vertex, position));
        glVertexArrayAttribFormat(vao, 1, 3, GL_FLOAT, GL_FALSE, offsetof(Vertex, normal));

        glVertexArrayAttribBinding(vao, 0, 0);
        glVertexArrayAttribBinding(vao, 1, 0);
    }

    ~DrawTrianglesShader()
    {
        if (vbo != 0) {
            glDeleteBuffers(1, &vbo);
        }
        if (vao != 0) {
            glDeleteVertexArrays(1, &vao);
        }
    }

    void Draw(
        const std::vector<Vertex>& vertices,
        const glm::mat4& MVP,
        const glm::mat4& Model,
        const glm::vec4& color,
        int objectId,
        const glm::vec3& lightDir = glm::normalize(glm::vec3(0.1f, 0.1f, -1.0f))
    )
    {
        if (vertices.empty())
            return;

        glNamedBufferData(vbo, static_cast<GLsizeiptr>(vertices.size() * sizeof(Vertex)), vertices.data(), GL_DYNAMIC_DRAW);

        use();
        SetUniformMat4("MVP", MVP);
        SetUniformMat4("Model", Model);
        SetUniform("Color", color);
        SetUniform("LightDir", lightDir);
        SetUniformI("ObjectId", objectId);

        glBindVertexArray(vao);
        glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(vertices.size()));
        glBindVertexArray(0);
        glUseProgram(0);
    }
};




class DrawAtomsShader : public Shader {
    static constexpr const char* vertexShaderSource = R"(
#version 430 core

struct RenderAtom {
    vec4 position; // {posX, posY, posZ, radius}
    vec4 color;    // {r, g, b, a}
    uvec4 flags;   // {x=highlight}
};

layout(std430, binding = 0) buffer RenderAtoms {
    RenderAtom atoms[];
};

uniform mat4 View;
uniform mat4 Proj;
uniform int  numVerticesPerAtom;
uniform float pi = 3.14159265359f;

out vec4 vertexColor;
flat out int atomId;
flat out uint highlight;
out vec2 localCoord;

void main() {
    const int numTrianglesPerAtom = numVerticesPerAtom - 2;
    float angle = 2.0f * pi * float(gl_VertexID) / float(numTrianglesPerAtom);

    vec4 atomPos = atoms[gl_InstanceID].position;
    atomId = int(atoms[gl_InstanceID].flags.y);

    vec4 viewSpacePos = View * vec4(atomPos.xyz, 1.0);
    float radius = atomPos.w;

    vec3 viewDir = normalize(-viewSpacePos.xyz);
    vec3 up = (abs(viewDir.z) < 0.999f) ? vec3(0.0, 0.0, 1.0) : vec3(0.0, 1.0, 0.0);
    vec3 right = normalize(cross(up, viewDir));
    vec3 up2   = cross(viewDir, right);

    vec4 posVS;
    float light;

    if (gl_VertexID == 0) {
        posVS = viewSpacePos;
        light = 0.7f;
        localCoord = vec2(0.0, 0.0);
    } else {
        float coneSlope = 0.15f;
        float coneDepth = radius * coneSlope;

        vec2 circle = vec2(cos(angle), sin(angle));
        vec3 offset3 =
            right * (circle.x * radius) +
            up2   * (circle.y * radius) -
            viewDir * coneDepth;

        posVS = viewSpacePos + vec4(offset3, 0.0);

        float ny = clamp(offset3.y / radius, -1.0f, 1.0f);
        light = clamp(ny * 0.5f + 0.6f, 0.0f, 1.0f);
        localCoord = circle;
    }

    highlight = atoms[gl_InstanceID].flags.x;
    vertexColor = vec4(atoms[gl_InstanceID].color.xyz * light, atoms[gl_InstanceID].color.w);
    gl_Position = Proj * posVS;
}
)";

    static constexpr const char* fragmentShaderSource = R"(
#version 430 core

in vec4 vertexColor;
flat in int atomId;
flat in uint highlight;
in vec2 localCoord;

layout(location = 0) out vec4 FragColor;
layout(location = 1) out int  FragAtomId;

void main() {
    vec3 color = vertexColor.rgb;

    if (highlight == 1u) {

        const float thickness = 0.08;
        const float intensity = 5.f;

        float r = length(localCoord);

        // thin bright ring near the outer edge        
        float halo = smoothstep(1.f-thickness, 1.f, r);
        
        // optional sharper falloff so it stays a rim instead of a wash
        halo *= halo * halo;

        if (halo > 0)
            color += color * vec3(halo*intensity);
    }

    FragColor = vec4(min(color, vec3(1.0f)), vertexColor.a);
    FragAtomId = atomId;
}
)";


    static constexpr int numVerticesPerAtom = 24;

public:
    DrawAtomsShader()
        : Shader(vertexShaderSource, fragmentShaderSource)
    {
    }

    ~DrawAtomsShader() {
    }

    void Draw(const SSBO& renderAtomsBuffer, int nAtoms, const glm::mat4& view, const glm::mat4& projection) {
        use();
        renderAtomsBuffer.Bind(0);

        SetUniformMat4("View", view);
        SetUniformMat4("Proj", projection);
        SetUniformI("numVerticesPerAtom", numVerticesPerAtom);

        glDrawArraysInstanced(GL_TRIANGLE_FAN, 0, numVerticesPerAtom, nAtoms);
        glUseProgram(0);
    }
};

class DrawAtomsPrettyShader : public Shader {
    struct SphereVertex {
        glm::vec3 position;
        glm::vec3 normal;
    };

    struct Triangle {
        uint32_t a;
        uint32_t b;
        uint32_t c;
    };

    GLuint vao = 0;
    GLuint vbo = 0;
    GLuint ebo = 0;
    GLsizei indexCount = 0;

    static constexpr const char* vertexShaderSource = R"(
#version 430 core

struct RenderAtom {
    vec4 position; // {posX_nm, posY_nm, posZ_nm, radius_nm}
    vec4 color;    // {r, g, b, a}
    uvec4 flags;   // {x=highlight}
};

layout(std430, binding = 0) buffer RenderAtoms {
    RenderAtom atoms[];
};

layout(location = 0) in vec3 inPosition;
layout(location = 1) in vec3 inNormal;

uniform mat4 View;
uniform mat4 Proj;

out vec3 fragNormalView;
out vec3 fragPositionView;
out vec4 vertexColor;
flat out int atomId;
flat out uint highlight;

void main() {
    RenderAtom atom = atoms[gl_InstanceID];

    vec3 worldPos = atom.position.xyz + inPosition * atom.position.w; // nm
    vec4 viewPos4 = View * vec4(worldPos, 1.0);

    fragPositionView = viewPos4.xyz;
    fragNormalView = normalize(mat3(View) * inNormal);

    vertexColor = atom.color;
    atomId = int(atoms[gl_InstanceID].flags.y);
    highlight = atom.flags.x;

    gl_Position = Proj * viewPos4;
}
)";

    static constexpr const char* fragmentShaderSource = R"(
#version 430 core

in vec3 fragNormalView;
in vec3 fragPositionView;
in vec4 vertexColor;
flat in int atomId;
flat in uint highlight;

layout(location = 0) out vec4 fragColor;
layout(location = 1) out int fragId;

float Luminance(vec3 color) {
    return dot(color, vec3(0.2126, 0.7152, 0.0722));
}

void main() {
    vec3 N = normalize(fragNormalView);
    vec3 V = normalize(-fragPositionView);

    // Light from above and slightly toward the camera in view space.
    vec3 L = normalize(vec3(0.18, 0.82, 0.55));
    vec3 H = normalize(L + V);

    vec3 baseColor = vertexColor.rgb;

    // Slight desaturation for a more muted, publication-like look.
    baseColor = mix(baseColor, vec3(Luminance(baseColor)), 0.10);

    float diffuse = max(dot(N, L), 0.0);
    float topBias = smoothstep(-0.15, 0.95, N.y);
    float bottomShade = smoothstep(0.05, 0.95, -N.y);
    float fresnel = pow(1.0 - max(dot(N, V), 0.0), 2.2);
    float specular = pow(max(dot(N, H), 0.0), 28.0);

    float lighting =
        0.34
        + 0.30 * diffuse
        + 0.34 * topBias;

    lighting *= (1.0 - 0.18 * bottomShade);

    vec3 color = baseColor * lighting;

    // Soft, broad specular. Keep it subtle.
    color += vec3(1.0) * specular * 0.10;

    // Soft rim/fresnel.
    color += vec3(1.0) * fresnel * 0.08;

    if (highlight != 0u) {
        vec3 highlightColor = vec3(1.0, 0.82, 0.32);
        color += highlightColor * fresnel * 0.45;
        color += highlightColor * specular * 0.25;
    }

    fragColor = vec4(color, vertexColor.a);
    fragId = atomId;
}
)";

public:
    DrawAtomsPrettyShader()
        : Shader(vertexShaderSource, fragmentShaderSource)
    {
        _CreateSphereMesh();
    }

    ~DrawAtomsPrettyShader() {
        if (ebo)
            glDeleteBuffers(1, &ebo);
        if (vbo)
            glDeleteBuffers(1, &vbo);
        if (vao)
            glDeleteVertexArrays(1, &vao);
    }

    void Draw(const SSBO& renderAtomsBuffer, int nAtoms, const glm::mat4& View, const glm::mat4& Proj) {
        if (nAtoms <= 0)
            return;

        use();

        renderAtomsBuffer.Bind(0);

        SetUniformMat4("View", View);
        SetUniformMat4("Proj", Proj);

        glBindVertexArray(vao);
        glDrawElementsInstanced(GL_TRIANGLES, indexCount, GL_UNSIGNED_INT, nullptr, nAtoms);

        glBindVertexArray(0);
    }

private:
    static uint64_t _MakeEdgeKey(uint32_t a, uint32_t b) {
        const uint32_t lo = std::min(a, b);
        const uint32_t hi = std::max(a, b);
        return (static_cast<uint64_t>(lo) << 32) | static_cast<uint64_t>(hi);
    }

    static uint32_t _GetMidpointIndex(
        uint32_t a,
        uint32_t b,
        std::vector<glm::vec3>& positions,
        std::unordered_map<uint64_t, uint32_t>& midpointCache
    ) {
        const uint64_t key = _MakeEdgeKey(a, b);

        auto it = midpointCache.find(key);
        if (it != midpointCache.end())
            return it->second;

        const glm::vec3 midpoint = glm::normalize((positions[a] + positions[b]) * 0.5f);
        const uint32_t index = static_cast<uint32_t>(positions.size());

        positions.push_back(midpoint);
        midpointCache.emplace(key, index);

        return index;
    }

    void _CreateSphereMesh() {
        std::vector<SphereVertex> vertices;
        std::vector<uint32_t> indices;
        _GenerateIcosphere(vertices, indices, 2); // Increase to 3 if you want it prettier/heavier.

        indexCount = static_cast<GLsizei>(indices.size());

        glGenVertexArrays(1, &vao);
        glGenBuffers(1, &vbo);
        glGenBuffers(1, &ebo);

        glBindVertexArray(vao);

        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glBufferData(
            GL_ARRAY_BUFFER,
            static_cast<GLsizeiptr>(vertices.size() * sizeof(SphereVertex)),
            vertices.data(),
            GL_STATIC_DRAW
        );

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
        glBufferData(
            GL_ELEMENT_ARRAY_BUFFER,
            static_cast<GLsizeiptr>(indices.size() * sizeof(uint32_t)),
            indices.data(),
            GL_STATIC_DRAW
        );

        glEnableVertexAttribArray(0);
        glVertexAttribPointer(
            0,
            3,
            GL_FLOAT,
            GL_FALSE,
            sizeof(SphereVertex),
            reinterpret_cast<void*>(offsetof(SphereVertex, position))
        );

        glEnableVertexAttribArray(1);
        glVertexAttribPointer(
            1,
            3,
            GL_FLOAT,
            GL_FALSE,
            sizeof(SphereVertex),
            reinterpret_cast<void*>(offsetof(SphereVertex, normal))
        );

        glBindVertexArray(0);
    }

    static void _GenerateIcosphere(
        std::vector<SphereVertex>& outVertices,
        std::vector<uint32_t>& outIndices,
        int subdivisions
    ) {
        const float t = (1.0f + std::sqrt(5.0f)) * 0.5f;

        std::vector<glm::vec3> positions = {
            glm::normalize(glm::vec3(-1.0f,  t, 0.0f)),
            glm::normalize(glm::vec3(1.0f,  t, 0.0f)),
            glm::normalize(glm::vec3(-1.0f, -t, 0.0f)),
            glm::normalize(glm::vec3(1.0f, -t, 0.0f)),

            glm::normalize(glm::vec3(0.0f, -1.0f,  t)),
            glm::normalize(glm::vec3(0.0f,  1.0f,  t)),
            glm::normalize(glm::vec3(0.0f, -1.0f, -t)),
            glm::normalize(glm::vec3(0.0f,  1.0f, -t)),

            glm::normalize(glm::vec3(t, 0.0f, -1.0f)),
            glm::normalize(glm::vec3(t, 0.0f,  1.0f)),
            glm::normalize(glm::vec3(-t, 0.0f, -1.0f)),
            glm::normalize(glm::vec3(-t, 0.0f,  1.0f))
        };

        std::vector<Triangle> triangles = {
            {0, 11, 5}, {0, 5, 1},  {0, 1, 7},  {0, 7, 10}, {0, 10, 11},
            {1, 5, 9},  {5, 11, 4}, {11, 10, 2},{10, 7, 6}, {7, 1, 8},
            {3, 9, 4},  {3, 4, 2},  {3, 2, 6},  {3, 6, 8},  {3, 8, 9},
            {4, 9, 5},  {2, 4, 11}, {6, 2, 10}, {8, 6, 7},  {9, 8, 1}
        };

        for (int i = 0; i < subdivisions; ++i) {
            std::unordered_map<uint64_t, uint32_t> midpointCache;
            std::vector<Triangle> nextTriangles;
            nextTriangles.reserve(triangles.size() * 4);

            for (const Triangle& tri : triangles) {
                const uint32_t ab = _GetMidpointIndex(tri.a, tri.b, positions, midpointCache);
                const uint32_t bc = _GetMidpointIndex(tri.b, tri.c, positions, midpointCache);
                const uint32_t ca = _GetMidpointIndex(tri.c, tri.a, positions, midpointCache);

                nextTriangles.push_back({ tri.a, ab, ca });
                nextTriangles.push_back({ tri.b, bc, ab });
                nextTriangles.push_back({ tri.c, ca, bc });
                nextTriangles.push_back({ ab, bc, ca });
            }

            triangles = std::move(nextTriangles);
        }

        outVertices.clear();
        outVertices.reserve(positions.size());

        for (const glm::vec3& p : positions) {
            SphereVertex v;
            v.position = p; // unit sphere
            v.normal = p;   // same on a unit sphere
            outVertices.push_back(v);
        }

        outIndices.clear();
        outIndices.reserve(triangles.size() * 3);

        for (const Triangle& tri : triangles) {
            outIndices.push_back(tri.a);
            outIndices.push_back(tri.b);
            outIndices.push_back(tri.c);
        }
    }
};
