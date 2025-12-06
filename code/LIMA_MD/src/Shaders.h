// No pragma, once include in Display.cpp
#include <GL/glew.h>

#include "MoleculeHull.cuh"
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

class DrawTrianglesShader : public Shader {
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


    DrawTrianglesShader() : Shader(hullVertexShaderSource, hullFragmentShaderSource) {
        // Generate and bind VAO
        glGenVertexArrays(1, &VAO);
        glBindVertexArray(VAO);
        glBindVertexArray(0);  // Unbind VAO
    }

    ~DrawTrianglesShader() {
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

template <bool isCUDA>
class DrawAtomsShader : public Shader {
    static constexpr const char* vertexShaderSource = R"(
#version 430 core 
struct RenderAtom { 
    vec4 position; // {posX, posY, posZ, radius} 
    vec4 color;    // {r, g, b, a} 
}; 
 
layout(std430, binding = 0) buffer RenderAtoms { 
    RenderAtom atoms[]; 
}; 

uniform mat4 View;   // View matrix
uniform mat4 Proj;   // Projection matrix
uniform int numAtoms; 
uniform int numVerticesPerAtom; 
uniform float pi = 3.14159265359f; 

out vec4 vertexColor; 

void main() { 
    int numTrianglesPerAtom = numVerticesPerAtom - 2; 
    float angle = 2.0f * pi * float(gl_VertexID) / float(numTrianglesPerAtom); 
    vec4 atomPos = atoms[gl_InstanceID].position; 

    // Atom center in view space
    vec4 viewSpacePos = View * vec4(atomPos.xyz, 1.0); 
    float radius = atomPos.w;

    // Build a basis so the cone axis points toward the camera (view-space)
    vec3 viewDir = normalize(-viewSpacePos.xyz);
    vec3 up = (abs(viewDir.z) < 0.999f)
        ? vec3(0.0, 0.0, 1.0)
        : vec3(0.0, 1.0, 0.0);
    vec3 right = normalize(cross(up, viewDir));
    vec3 up2   = cross(viewDir, right);

    vec4 objectSpacePos;
    float light;

    if (gl_VertexID == 0) {
        // Center vertex at atom center
        objectSpacePos = viewSpacePos;
        light = .7f;
    } else {
        // Cone geometry in view space
        float coneSlope = 0.15f;
        float coneDepth = radius * coneSlope;

        vec3 offset3 =
            right * (cos(angle) * radius) +
            up2   * (sin(angle) * radius) -
            viewDir * coneDepth;

        objectSpacePos = viewSpacePos + vec4(offset3, 0.0);
        

        // Use vertical (view-space Y) offset on the circle.
        float ny = clamp(offset3.y / radius, -1.0f, 1.0f); // -1 bottom, +1 top
        // Map ny from [-1,1] -> [0,1], then add some ambient floor
        //light = 0.3f + 0.7f * (0.5f + 0.5f * ny);
        float light1 = ((ny * 0.5f) * (ny * 0.5f) + 0.5f);
        light = clamp(ny * 0.5f + 0.6, 0.0f, 1.0f);
    }

    gl_Position = Proj * objectSpacePos;

    vertexColor = vec4(atoms[gl_InstanceID].color.xyz * light,
                       atoms[gl_InstanceID].color.w); 
}  
    )";


    static constexpr const char* fragmentShaderSource = R"(
#version 430 core 
in vec4 vertexColor; 
out vec4 FragColor; 
  
void main() { 
    FragColor = vertexColor; 
} 


//in vec4 vertexColor;
//in vec2 localUV;
//
//out vec4 FragColor;
//
//void main() {
//    float r2 = dot(localUV, localUV);
//    if (r2 > 1.0)
//        discard;
//
//    // simple radial attenuation (edge darker)
//    float radial = 1.0 - r2;
//
//    FragColor = vec4(vertexColor.rgb * radial, vertexColor.a);
//}



//    #version 430 core 
//    in vec4 vertexColor; 
//    out vec4 FragColor; 
//    in vec2 localUV;
//
//    void main() { 
////        FragColor = vertexColor; 
//        float r2 = dot(localUV, localUV);
//        if (r2 > 1.0) discard;           // HARD sphere silhouette
//
//        float z = sqrt(1.0 - r2);        // sphere depth
//        vec3 normal = normalize(vec3(localUV, z));
//
//        float light = max(normal.z, 0.0);
//        FragColor = vec4(vertexColor.rgb * light, vertexColor.a);
//
//        /*
//        vec3 lightDir = normalize(vec3(0.3, 0.5, 0.8)); // view-space
//        float NdotL = max(dot(normal, lightDir), 0.0);
//
//        float ambient = 0.2;
//        vec3 color = vertexColor.rgb * (ambient + 0.8 * NdotL);
//
//        FragColor = vec4(color, vertexColor.a);*/        
//    } 
    )";

    GLuint VBO;
public:
    SSBO renderAtomsBuffer{};

    const int numAtomsReservedInRenderatomsBuffer;

    DrawAtomsShader(int numAtoms, cudaGraphicsResource** renderAtomsBufferSource) : 
        Shader(vertexShaderSource, fragmentShaderSource),
        numAtomsReservedInRenderatomsBuffer(numAtoms)
    {

    // Generate and bind VBO
    glGenBuffers(1, &VBO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);

    if constexpr (isCUDA) {
        // Allocate the renderAtomsbuffer, and register it with CUDA
        renderAtomsBuffer.Resize(numAtoms * sizeof(RenderAtom));
        cudaGraphicsGLRegisterBuffer(renderAtomsBufferSource, renderAtomsBuffer.GetID(), cudaGraphicsMapFlagsWriteDiscard);
    }
    //Otherwise, whoever preps will write directly to our renderAtomsBuffer. That is bad tho, they should maybe just call another constructor of my class with the SSBO ready?
    // Allocate space for the VBO. We assume `sizeof(RenderAtom)` is the size of the data.
    glBufferData(GL_ARRAY_BUFFER, numAtoms * sizeof(RenderAtom), nullptr, GL_DYNAMIC_DRAW);

    // Enable and set up position attribute
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, sizeof(RenderAtom), (void*)offsetof(RenderAtom, position));

    // Enable and set up color attribute
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(RenderAtom), (void*)offsetof(RenderAtom, color));

    // Unbind VAO and VBO to prevent unintended modifications
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    ~DrawAtomsShader() {
		glDeleteBuffers(1, &VBO);
	}
    // This functions needs no renderAtoms input, as the class is initialized with a reference to the buffer where
    // CUDA will place the input data
    void Draw(const glm::mat4& view, const glm::mat4& projection, int nAtoms) {

        // If we call this shader multiple times it doesnt matter if the following times has less atoms, but we currently have
        // not implemented to ability to increase the buffer size
        if (nAtoms > numAtomsReservedInRenderatomsBuffer) {
			throw std::runtime_error("Number of atoms in DrawAtomsShader::Draw does not match the number of atoms in the constructor. \n\
                This is likely due to the Renderer being tasked with both rendering compounds and moleculeHulls in the same instance");
		}

        use();

        renderAtomsBuffer.Bind(0);

        //SetUniformMat4("MVP", MVP);
        SetUniformMat4("View", view);
        SetUniformMat4("Proj", projection);
        SetUniformI("numAtoms", nAtoms);

        const int numVerticesPerAtom = 24;
		SetUniformI("numVerticesPerAtom", numVerticesPerAtom);
        
        glDrawArraysInstanced(GL_TRIANGLE_FAN, 0, numVerticesPerAtom, nAtoms);

        glBindVertexArray(0);
        glUseProgram(0);
    }

};