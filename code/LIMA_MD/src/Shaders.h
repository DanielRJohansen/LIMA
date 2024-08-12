// No pragma, once include in Display.cpp
#include <GL/glew.h>
#include "GLShader.h"
#include <glm.hpp>
#include <gtc/matrix_transform.hpp>
#include <gtc/type_ptr.hpp>


#include <cuda_gl_interop.h>


class DrawBoxOutlineShader : public Shader {
    static constexpr const char* vertexShaderSource = R"(
        #version 430 core
        layout(location = 0) in vec3 aPos;
        layout(location = 1) in vec3 aColor;
        out vec3 vertexColor;
        uniform mat4 MVP;
        void main() {
            gl_Position = MVP * vec4(aPos, 1.0);
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
            -0.5f, -0.5f, -0.5f,  0.2f, 0.2f, 0.8f,
             0.5f, -0.5f, -0.5f,  0.2f, 0.2f, 0.8f,
             0.5f,  0.5f, -0.5f,  0.2f, 0.2f, 0.8f,
            -0.5f,  0.5f, -0.5f,  0.2f, 0.2f, 0.8f,
            -0.5f, -0.5f,  0.5f,  0.4f, 0.4f, 0.4f,
             0.5f, -0.5f,  0.5f,  0.4f, 0.4f, 0.4f,
             0.5f,  0.5f,  0.5f,  0.8f, 0.4f, 0.4f,
            -0.5f,  0.5f,  0.5f,  0.8f, 0.4f, 0.4f,
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

    void Draw(const glm::mat4 MVP) {
        use();

        glBindVertexArray(VAO);
        SetUniformMat4("MVP", MVP);
        glDrawElements(GL_LINES, 24, GL_UNSIGNED_INT, 0);
        glBindVertexArray(0);

        glUseProgram(0);
    }
};

class DrawTrianglesShader : public Shader {
    static constexpr const char* hullVertexShaderSource = R"(
    #version 430 core

    struct Tri {
        float data[12]; // DONT CHANGE, it's OpenGL being annoying with offsets
        // vec3[3] vertices;
        // vec3 normal
    };

    layout(std430, binding = 0) buffer TriBuffer {
        Tri tris[];
    };

    uniform mat4 MVP;
    out vec3 fragColor;

    uniform vec3 lightDir = vec3(0.0, 0.0, -1.0); // Light coming from directly above
    
    uniform int drawMode; // 0 = FACES, 1 = EDGES
    

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
        vec3 position = vec3(tris[triIndex].data[vertexIndex * 3], 
                             tris[triIndex].data[vertexIndex * 3 + 1], 
                             tris[triIndex].data[vertexIndex * 3 + 2]);
        gl_Position = MVP * vec4(position, 1.0);



        // Now set color
        vec3 triNormal = vec3(tris[triIndex].data[3*3+0], 
                        tris[triIndex].data[3*3+1], 
                        tris[triIndex].data[3*3+2]);

        // Simulate area light by blending the normal with the light direction
        float brightness = clamp(
            dot(-triNormal, lightDir) * 0.5f + 0.5f, 
            0.1f,
            1.f);
        
        if (drawMode == 1) { // For EDGES mode, color the edges red  			
			fragColor = vec3(1.0, 0.0, 0.0);
		}
        else { // Generate a pseudo-random color based on the triangle index        
            fragColor = brightness * GenerateRandomColor(triIndex);
        }
    }
)";


    static constexpr const char* hullFragmentShaderSource = R"(
    #version 430 core

    in vec3 fragColor;
    out vec4 color;

    void main() {
        color = vec4(fragColor, 1.0);
    }
)";

    GLuint VAO, SSBO;

public:
    enum DrawMode {
        FACES=0,
        EDGES=1
    };

    DrawTrianglesShader() : Shader(hullVertexShaderSource, hullFragmentShaderSource) {
        // Generate and bind VAO
        glGenVertexArrays(1, &VAO);
        glBindVertexArray(VAO);

        // Generate SSBO (Shader Storage Buffer Object)
        glGenBuffers(1, &SSBO);

        glBindVertexArray(0);  // Unbind VAO
    }

    ~DrawTrianglesShader() {
        glDeleteVertexArrays(1, &VAO);
        glDeleteBuffers(1, &SSBO);
    }

    void Draw(const glm::mat4& MVP, const std::vector<Tri>& tris, DrawMode mode) {
        use(); // Use the shader program

        // Update SSBO with new triangle data
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, SSBO);
        glBufferData(GL_SHADER_STORAGE_BUFFER, tris.size() * sizeof(Tri), tris.data(), GL_DYNAMIC_DRAW);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, SSBO);  // Binding index 0 matches the shader
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

        // Set the MVP matrix uniform
        SetUniformMat4("MVP", MVP);

        // Set the draw mode uniform
        SetUniformI("drawMode", mode == FACES ? 0 : 1);

        // Draw the triangles based on the selected mode
        glBindVertexArray(VAO);
        if (mode == FACES) {
            glDrawArrays(GL_TRIANGLES, 0, tris.size() * 3);
        }
        else if (mode == EDGES) {
            glDrawArrays(GL_LINES, 0, tris.size() * 6); // Each triangle has 3 edges, each edge has 2 vertices
        }
        glBindVertexArray(0);
        glUseProgram(0);  // Unbind shader program
    }
};




class DrawAtomsShader : public Shader {
    static constexpr const char* vertexShaderSource = R"(
    #version 430 core 
    struct RenderAtom { 
    vec4 position; // {posX, posY, posZ, radius} 
    vec4 color;    // {r, g, b, a} 
    }; 
     
    // Uniform block for RenderAtoms, dynamic size 
    layout(std430, binding = 0) buffer RenderAtoms { 
        RenderAtom atoms[]; 
    }; 
  
    uniform mat4 MVP; 
    uniform int numAtoms; 
    uniform int numVerticesPerAtom; 
    uniform float pi = 3.14159265359f; 
     
    out vec4 vertexColor; 
     
    void main() { 
        int numTrianglesPerAtom = numVerticesPerAtom - 2; 
        float light = (sin(float(gl_VertexID * 2 * pi) / float(numTrianglesPerAtom )) + 1.f) / 2.f;
        float angle = 2.0f * pi * float(gl_VertexID) / float(numTrianglesPerAtom); 
        vec4 atomPos = atoms[gl_InstanceID].position; 
        vec4 worldPos = vec4(atomPos.xyz, 1.0); 
        vec4 offset = vec4(cos(angle) * atomPos.w, sin(angle)  * atomPos.w, 0.5, 0.0) * 0.3f; 
        gl_Position = MVP * worldPos + offset; 
        if (gl_VertexID == 0) { 
            gl_Position = MVP * vec4(atomPos.xyz, 1.0); 
        } 
        vertexColor = atoms[gl_InstanceID].color * light; 
    } 
    )";


    static constexpr const char* fragmentShaderSource = R"(
    #version 430 core 
    in vec4 vertexColor; 
    out vec4 FragColor; 
     
    void main() { 
        FragColor = vertexColor; 
    } 
    )";

    GLuint VBO;

    public
        DrawAtomsShader(int numAtoms)
};