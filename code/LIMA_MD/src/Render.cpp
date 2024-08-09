#include "DisplayV2.h"

#include <glm.hpp>
#include <gtc/matrix_transform.hpp>
#include <gtc/type_ptr.hpp>


#include <cuda_gl_interop.h>



// ------------------------------------------------------------- Render functions used during sim ------------------------------------------------------------- //

glm::mat4 GetMVPMatrix(float camera_distance, float camera_pitch, float camera_yaw, int screenWidth, int screenHeight) {
    // Set up the perspective projection
    double aspectRatio = static_cast<double>(screenWidth) / static_cast<double>(screenHeight);
    double fovY = 45.0;
    double nearPlane = 0.1;
    double farPlane = 10.0;
    double fH = tan(fovY / 360.0 * 3.141592653589793) * nearPlane;
    double fW = fH * aspectRatio;

    glm::mat4 projection = glm::frustum(-fW, fW, -fH, fH, nearPlane, farPlane);

    // Set up the view matrix
    glm::mat4 view = glm::mat4(1.0f);
    view = glm::translate(view, glm::vec3(0.0f, 0.0f, camera_distance));
    view = glm::rotate(view, glm::radians(-90.0f), glm::vec3(1.0f, 0.0f, 0.0f));  // Fixed rotation around x-axis
    view = glm::rotate(view, glm::radians(camera_pitch), glm::vec3(1.0f, 0.0f, 0.0f));  // Rotation around x-axis for pitch
    view = glm::rotate(view, glm::radians(camera_yaw), glm::vec3(0.0f, 0.0f, 1.0f));  // Rotation around z-axis for yaw (y-axis in GLM's coordinate system)

    // Model matrix (identity in this case)
    glm::mat4 model = glm::mat4(1.0f);

    return projection * view * model;
}




void checkCompileErrors(GLuint shader, std::string type) {
    GLint success;
    GLchar infoLog[1024];
    if (type != "PROGRAM") {
        glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
        if (!success) {
            glGetShaderInfoLog(shader, 1024, NULL, infoLog);
            std::cerr << "\n| ERROR::SHADER-COMPILATION-ERROR of type: " << type << "\n" << infoLog << "\n -- --------------------------------------------------- -- " << std::endl;
        }
    }
    else {
        glGetProgramiv(shader, GL_LINK_STATUS, &success);
        if (!success) {
            glGetProgramInfoLog(shader, 1024, NULL, infoLog);
            std::cerr << "\n| ERROR::PROGRAM-LINKING-ERROR of type: " << type << "\n" << infoLog << "\n -- --------------------------------------------------- -- " << std::endl;
        }
    }
}




//GLuint CompileShader(GLenum type, const char* source) {
//    GLuint shader = glCreateShader(type);
//    glShaderSource(shader, 1, &source, nullptr);
//    glCompileShader(shader);
//
//    // Error checking
//    GLint success;
//    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
//    if (!success) {
//        // Retrieve and print the error message
//        char infoLog[512];
//        glGetShaderInfoLog(shader, 512, nullptr, infoLog);
//        std::string shaderType = (type == GL_VERTEX_SHADER) ? "VERTEX" : "FRAGMENT";
//        throw std::runtime_error(shaderType + " SHADER COMPILATION FAILED:\n" + infoLog);
//    }
//    return shader;
//}
//
//GLuint CreateShaderProgram() {
//    GLuint vertexShader = CompileShader(GL_VERTEX_SHADER, vertexShaderSource);
//    GLuint fragmentShader = CompileShader(GL_FRAGMENT_SHADER, fragmentShaderSource);
//
//    // Link shaders into a program
//    GLuint shaderProgram = glCreateProgram();
//    glAttachShader(shaderProgram, vertexShader);
//    glAttachShader(shaderProgram, fragmentShader);
//    glLinkProgram(shaderProgram);
//
//    // Error checking
//    GLint success;
//    glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
//    if (!success) {
//        // Retrieve and print the error message
//        char infoLog[512];
//        glGetProgramInfoLog(shaderProgram, 512, nullptr, infoLog);
//        throw std::runtime_error("SHADER PROGRAM LINKING FAILED:\n" + infoLog);
//    }
//
//    // Cleanup shaders as they're linked now
//    glDeleteShader(vertexShader);
//    glDeleteShader(fragmentShader);
//
//    return shaderProgram;
//}
//















const char* box_vertexShaderSource = " \
#version 430 core \n\
layout(location = 0) in vec3 aPos; \n\
layout(location = 1) in vec3 aColor; \n\
 \n\
out vec3 vertexColor; \n\
 \n\
uniform mat4 MVP; // Model,View,Perspective \n\
 \n\
void main() { \n\
    gl_Position = MVP * vec4(aPos, 1.0); \n\
    vertexColor = aColor; \n\
} \n\
";

const char* box_fragmentShaderSource = " \
#version 430 core \n\
in vec3 vertexColor; \n\
out vec4 FragColor; \n\
 \n\
void main() { \n\
    FragColor = vec4(vertexColor, 1.0); \n\
} \n\
";

GLuint boxVAO, boxVBO, boxEBO;

void SetupDrawBoxPipeline(GLuint& drawBoxShaderProgram) {
    // Compile and link box shaders
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &box_vertexShaderSource, nullptr);
    glCompileShader(vertexShader);
    checkCompileErrors(vertexShader, "VERTEX");
    

    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &box_fragmentShaderSource, nullptr);
    glCompileShader(fragmentShader);
    checkCompileErrors(fragmentShader, "FRAGMENT");

    drawBoxShaderProgram = glCreateProgram();
    glAttachShader(drawBoxShaderProgram, vertexShader);
    glAttachShader(drawBoxShaderProgram, fragmentShader);
    glLinkProgram(drawBoxShaderProgram);
    checkCompileErrors(drawBoxShaderProgram, "PROGRAM");

    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

    // Box vertices and colors
    float boxVertices[] = {
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

    unsigned int boxIndices[] = {
        // bottom
        0, 1,
        1, 2,
        2, 3,
        3, 0,
        // top
        4, 5,
        5, 6,
        6, 7,
        7, 4,
        // sides
        0, 4,
        1, 5,
        2, 6,
        3, 7
    };

    const GLsizei numVAOs = 1;
    const GLsizei numVBOs = 1;
    const GLsizei numEBOs = 1;

    const GLint positionAttributeIndex = 0;
    const GLint colorAttributeIndex = 1;
    const GLint positionComponentCount = 3;
    const GLint colorComponentCount = 3;
    const GLsizei stride = (positionComponentCount + colorComponentCount) * sizeof(float);
    const void* positionOffset = (void*)0;
    const void* colorOffset = (void*)(positionComponentCount * sizeof(float));

    glGenVertexArrays(numVAOs, &boxVAO);
    glGenBuffers(numVBOs, &boxVBO);
    glGenBuffers(numEBOs, &boxEBO);

    glBindVertexArray(boxVAO);

    glBindBuffer(GL_ARRAY_BUFFER, boxVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(boxVertices), boxVertices, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, boxEBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(boxIndices), boxIndices, GL_STATIC_DRAW);

    glVertexAttribPointer(positionAttributeIndex, positionComponentCount, GL_FLOAT, GL_FALSE, stride, positionOffset);
    glEnableVertexAttribArray(positionAttributeIndex);
    glVertexAttribPointer(colorAttributeIndex, colorComponentCount, GL_FLOAT, GL_FALSE, stride, colorOffset);
    glEnableVertexAttribArray(colorAttributeIndex);

    glBindVertexArray(0);
}


void Display::DrawBoxOutline() {
    glUseProgram(*drawBoxShaderProgram);

    const GLuint mvpLocation = glGetUniformLocation(*drawBoxShaderProgram, "MVP");
    const glm::mat4 MVP = GetMVPMatrix(camera_distance, camera_pitch * rad2deg, camera_yaw * rad2deg, screenWidth, screenHeight);
    glUniformMatrix4fv(mvpLocation, 1, GL_FALSE, glm::value_ptr(MVP));

    const GLsizei boxIndexCount = 24;
    glBindVertexArray(boxVAO);
    glDrawElements(GL_LINES, boxIndexCount, GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);

    glUseProgram(0);
}

















const char* vertexShaderSource =
" \
#version 430 core \n\
struct RenderAtom { \n\
vec4 position; // {posX, posY, posZ, radius} \n\
vec4 color;    // {r, g, b, a} \n\
}; \n\
 \n\
// Uniform block for RenderAtoms, dynamic size \n\
layout(std430, binding = 0) buffer RenderAtoms { \n\
    RenderAtom atoms[]; \n\
}; \n\
\
uniform mat4 MVP; \n\
uniform int numAtoms; \n\
uniform int numVerticesPerAtom; \n\
uniform float pi = 3.14159265359f; \n\
 \n\
out vec4 vertexColor; \n\
 \n\
void main() { \n\
    int numTrianglesPerAtom = numVerticesPerAtom - 2; \n\
    float light = (sin(float(gl_VertexID * 2 * pi) / float(numTrianglesPerAtom )) + 1.f) / 2.f;\n\
    float angle = 2.0f * pi * float(gl_VertexID) / float(numTrianglesPerAtom); \n\
    vec4 atomPos = atoms[gl_InstanceID].position; \n\
    vec4 worldPos = vec4(atomPos.xyz, 1.0); \n\
    vec4 offset = vec4(cos(angle) * atomPos.w, sin(angle)  * atomPos.w, 0.5, 0.0) * 0.3f; \n\
    gl_Position = MVP * worldPos + offset; \n\
    if (gl_VertexID == 0) { \n\
        gl_Position = MVP * vec4(atomPos.xyz, 1.0); \n\
    } \n\
    vertexColor = atoms[gl_InstanceID].color * light; \n\
} \n\
";


const char* fragmentShaderSource = " \
#version 430 core \n\
in vec4 vertexColor; \n\
out vec4 FragColor; \n\
 \n\
void main() { \n\
    FragColor = vertexColor; \n\
} \n\
";








GLuint drawAtomsVBO;

void SetupDrawAtomsPipeline(size_t numAtoms, GLuint& drawAtomsShaderProgram, std::optional<GLuint>& renderAtomsBuffer, cudaGraphicsResource** renderAtomsBufferCUDA) {
    // Create and compile the vertex shader
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexShaderSource, nullptr);
    glCompileShader(vertexShader);
    checkCompileErrors(vertexShader, "VERTEX");

    // Create and compile the fragment shader
    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentShaderSource, nullptr);
    glCompileShader(fragmentShader);
    checkCompileErrors(fragmentShader, "FRAGMENT");

    // Link shaders into a shader program
    drawAtomsShaderProgram = glCreateProgram();
    glAttachShader(drawAtomsShaderProgram, vertexShader);
    glAttachShader(drawAtomsShaderProgram, fragmentShader);
    glLinkProgram(drawAtomsShaderProgram);
    checkCompileErrors(drawAtomsShaderProgram, "PROGRAM");

    // Clean up shaders as they are linked now
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

    // Generate and bind VBO
    glGenBuffers(1, &drawAtomsVBO);
    glBindBuffer(GL_ARRAY_BUFFER, drawAtomsVBO);

    // Generate SSBOs
    renderAtomsBuffer.emplace(0);
    glGenBuffers(1, &renderAtomsBuffer.value());
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, *renderAtomsBuffer);
    glBufferData(GL_SHADER_STORAGE_BUFFER, numAtoms * sizeof(RenderAtom), nullptr, GL_DYNAMIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
	cudaGraphicsGLRegisterBuffer(renderAtomsBufferCUDA, *renderAtomsBuffer, cudaGraphicsMapFlagsWriteDiscard);

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




void Display::DrawAtoms(size_t numAtoms) {
    const int numVerticesPerAtom = 12; // Number of vertices to define a disc (30 for the circle + 1 center + 1 to close the loop)

    glUseProgram(*drawAtomsShaderProgram);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, *renderAtomsBuffer);

    const GLuint mvpLocation = glGetUniformLocation(*drawAtomsShaderProgram, "MVP");
    const glm::mat4 MVP = GetMVPMatrix(camera_distance, camera_pitch * rad2deg, camera_yaw * rad2deg, screenWidth, screenHeight);
    glUniformMatrix4fv(mvpLocation, 1, GL_FALSE, glm::value_ptr(MVP));

    const GLuint numAtomsLocation = glGetUniformLocation(*drawAtomsShaderProgram, "numAtoms");
    glUniform1i(numAtomsLocation, numAtoms);

    const GLuint numVerticesPerAtomLocation = glGetUniformLocation(*drawAtomsShaderProgram, "numVerticesPerAtom");
    glUniform1i(numVerticesPerAtomLocation, numVerticesPerAtom);

    glDrawArraysInstanced(GL_TRIANGLE_FAN, 0, numVerticesPerAtom, numAtoms);

    glBindVertexArray(0);
    glUseProgram(0);
}



















const char* hullVertexShaderSource = R"(
#version 430 core
layout(location = 0) in vec4 inPosition;
//layout(location = 1) in vec3 inNormal;
//
//struct Tri {
//    vec3 vertices[3];
//    vec3 normal;
//};
//layout(std140, binding = 0) uniform TriBuffer {
//    Tri tris[];
//};

uniform mat4 MVP;

//out vec3 fragNormal;

void main() {
    // Assume you want to access the first triangle
    //Tri tri = tris[0];
    //
    //vec3 v0 = tri.vertices[0];
    //vec3 normal = tri.normal;
    //
    //// Process with the normal and vertices
    //gl_Position = MVP * vec4(v0, 1.0);


    gl_Position = MVP * vec4(inPosition.xyz, 1.0);
    //fragNormal = inNormal;
}
)";

const char* hullFragmentShaderSource = R"(
#version 430 core
//in vec3 fragNormal;
//out vec4 FragColor;
//uniform vec3 lightDir = normalize(vec3(1.0, 1.0, 1.0));
//uniform vec3 objectColor = vec3(0.4, 0.6, 0.8);

out vec4 color;

void main() {
    //float brightness = max(dot(normalize(fragNormal), lightDir), 0.0);
    //vec3 diffuse = brightness * objectColor;
    //FragColor = vec4(diffuse, 1.0);
	//FragColor = vec4(objectColor, 1.0);
    color = vec4(1,0,0,1);
}
)";

GLuint trisVAO, trisVBO, trisShaderProgram;

void SetupDrawTrisPipeline() {
    // Create and compile the vertex shader
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &hullVertexShaderSource, nullptr);
    glCompileShader(vertexShader);
    checkCompileErrors(vertexShader, "VERTEX");

    // Create and compile the fragment shader
    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &hullFragmentShaderSource, nullptr);
    glCompileShader(fragmentShader);
    checkCompileErrors(fragmentShader, "FRAGMENT");

    // Link shaders into a shader program
    trisShaderProgram = glCreateProgram();
    glAttachShader(trisShaderProgram, vertexShader);
    glAttachShader(trisShaderProgram, fragmentShader);
    glLinkProgram(trisShaderProgram);
    checkCompileErrors(trisShaderProgram, "PROGRAM");

    // Clean up shaders as they are linked now
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

    // Generate and bind VAO
    glGenVertexArrays(1, &trisVAO);
    glBindVertexArray(trisVAO);

    // Generate and bind VBO
    glGenBuffers(1, &trisVBO);
    glBindBuffer(GL_ARRAY_BUFFER, trisVBO);

    // Position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Float3), (void*)0);
    glEnableVertexAttribArray(0);

    //// Normal attribute
    //glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Float3), (void*)offsetof(Tri, normal));
    //glEnableVertexAttribArray(1);

    glBindVertexArray(0);  // Unbind VAO
}



















void Display::initializePipeline(size_t numAtoms) {
    drawBoxShaderProgram.emplace(0);
    drawAtomsShaderProgram.emplace(0);
    //drawMoleculeContainersProgram.emplace(0);

	// Set up the draw box pipeline
	SetupDrawBoxPipeline(*drawBoxShaderProgram);

	// Set up the draw atoms pipeline
	SetupDrawAtomsPipeline(numAtoms, *drawAtomsShaderProgram, renderAtomsBuffer, &renderAtomsBufferCudaResource);


    //SetupDrawTrisPipeline(*drawMoleculeContainersProgram);
    SetupDrawTrisPipeline();
    

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}






void Display::TerminateGLEW() {
    if (!drawAtomsShaderProgram.has_value())
        return; // Pipeline was never initialized, so nothing to clean up

    // Cleanup
    glDeleteVertexArrays(1, &drawAtomsVBO);
    glDeleteVertexArrays(1, &boxVAO);
    glDeleteVertexArrays(1, &boxVBO);
    glDeleteVertexArrays(1, &renderAtomsBuffer.value());
    

    if (drawAtomsShaderProgram.has_value())
        glDeleteProgram(*drawAtomsShaderProgram);
    if (drawBoxShaderProgram.has_value())
        glDeleteProgram(*drawBoxShaderProgram);
}








// ------------------------------ Render functions for various other tasks, less focus on performance ------------------------------ //

//void Display::DrawMoleculeContainers(const std::vector<MoleculeContainerSmall>& molecules, float boxlenNM) {
//    const glm::mat4 MVP = GetMVPMatrix(camera_distance, camera_pitch * rad2deg, camera_yaw * rad2deg, screenWidth, screenHeight);
//
//    //glBegin(GL_LINES);
//    //glColor3f(1.0f, 0, 0); // Set edge color to red
//
//    //for (const auto& molecule : molecules) {
//    //    const std::vector<Plane>& facets = molecule.convexHull.GetFacets();
//
//    //    for (const Plane& facet : facets) {
//    //        // Render edges of the triangular facet
//    //        for (size_t i = 0; i < 3; ++i) {
//    //            const Float3& v1 = facet.vertices[i] / boxlenNM - 0.5f;
//    //            const Float3& v2 = facet.vertices[(i + 1) % 3] / boxlenNM - 0.5f;
//
//    //            // Transform vertices using MVP matrix
//    //            glm::vec4 v1Transformed = MVP * glm::vec4(v1.x, v1.y, v1.z, 1.0f);
//    //            glm::vec4 v2Transformed = MVP * glm::vec4(v2.x, v2.y, v2.z, 1.0f);
//
//    //            glVertex3f(v1Transformed.x / v1Transformed.w, v1Transformed.y / v1Transformed.w, v1Transformed.z / v1Transformed.w);
//    //            glVertex3f(v2Transformed.x / v2Transformed.w, v2Transformed.y / v2Transformed.w, v2Transformed.z / v2Transformed.w);
//    //        }
//    //    }
//    //}
//    //glEnd();
//
//
//    glBegin(GL_TRIANGLES);
//    srand(42); // Set random seed to ensure consistent color
//    for (const auto& molecule : molecules) {
//        const auto& facets = molecule.convexHull.GetFacets();
//
//        for (const Plane& facet : facets) {
//            // Generate a random color with 10% transparency
//            glColor4f(static_cast<float>(rand()) / RAND_MAX,
//                static_cast<float>(rand()) / RAND_MAX,
//                static_cast<float>(rand()) / RAND_MAX,
//                1);
//
//            for (size_t i = 0; i < 3; ++i) {
//                const Float3& v = facet.vertices[i] / boxlenNM - 0.5f;
//
//                // Transform vertex using MVP matrix
//                glm::vec4 vTransformed = MVP * glm::vec4(v.x, v.y, v.z, 1.0f);
//
//                glVertex3f(vTransformed.x / vTransformed.w, vTransformed.y / vTransformed.w, vTransformed.z / vTransformed.w);
//            }
//        }
//    }
//    glEnd();
//
//
//
//
//
//    glBegin(GL_TRIANGLES);
//    for (const auto& molecule : molecules) {
//        for (const auto& pos : molecule.particlePositions) {
//            glm::vec4 atomPos = glm::vec4(pos.x / boxlenNM - 0.5f, pos.y / boxlenNM - 0.5f, pos.z / boxlenNM - 0.5f, 1.0f);
//            glm::vec4 transformedPos = MVP * atomPos;
//
//            // Draw the disc as a series of triangles
//            const int numVerticesPerAtom = 16;
//            const float pi = 3.14159265359f;
//
//            for (int i = 0; i < numVerticesPerAtom; ++i) {
//                float angle1 = 2.0f * pi * i / numVerticesPerAtom;
//                float angle2 = 2.0f * pi * (i + 1) / numVerticesPerAtom;
//
//                const float scale = 0.01f;
//
//                glm::vec4 offset1 = glm::vec4(cos(angle1) * scale, sin(angle1) * scale, 0.0f, 0.0f);
//                glm::vec4 offset2 = glm::vec4(cos(angle2) * scale, sin(angle2) * scale, 0.0f, 0.0f);
//
//                glm::vec4 v1 = transformedPos + offset1;
//                glm::vec4 v2 = transformedPos + offset2;
//
//
//                // Calculate light intensity based on angle
//                float light1 = (sin(angle1) + 1.0f) / 2.0f;
//                float light2 = (sin(angle2) + 1.0f) / 2.0f;
//
//                // Set color with varying intensity
//                glColor3f(1.0f * light1, 1.0f * light1, 0.0f * light1);
//                glVertex3f(transformedPos.x / transformedPos.w, transformedPos.y / transformedPos.w, transformedPos.z / transformedPos.w);
//
//                glColor3f(1.0f * light1, 1.0f * light1, 0.0f * light1);
//                glVertex3f(v1.x / v1.w, v1.y / v1.w, v1.z / v1.w);
//
//                glColor3f(1.0f * light2, 1.0f * light2, 0.0f * light2);
//                glVertex3f(v2.x / v2.w, v2.y / v2.w, v2.z / v2.w);
//            }
//        }
//    }
//    glEnd();
//
//}







void Display::DrawMoleculeContainers(const std::vector<MoleculeContainerSmall>& molecules, float boxlenNM) {
    const glm::mat4 MVP = GetMVPMatrix(camera_distance, camera_pitch * rad2deg, camera_yaw * rad2deg, screenWidth, screenHeight);

    for (const auto& molecule : molecules) {

        std::vector<Tri> tris(molecule.convexHull.numFacets);
        for (size_t i = 0; i < molecule.convexHull.numFacets; ++i) {
			const Plane& facet = molecule.convexHull.GetFacets()[i];
			Tri& tri = tris[i];

			for (size_t j = 0; j < 3; ++j) {
				const Float3& v = facet.vertices[j] / boxlenNM - 0.5f;
				tri.vertices[j] = v;
			}
//            tri.normal = facet.normal;
		}

        glUseProgram(trisShaderProgram);
        glBindVertexArray(trisVAO);
        glBindBuffer(GL_ARRAY_BUFFER, trisVBO);

        // Upload the tris data to the GPU
        glBufferData(GL_ARRAY_BUFFER, tris.size() * sizeof(Tri), tris.data(), GL_DYNAMIC_DRAW);

        // Set the MVP matrix uniform
        
        GLuint mvpLocation = glGetUniformLocation(trisShaderProgram, "MVP");
        glUniformMatrix4fv(mvpLocation, 1, GL_FALSE, glm::value_ptr(MVP));

        // Draw the triangles
        glDrawArrays(GL_TRIANGLES, 0, tris.size() * 3);

        glBindVertexArray(0);  // Unbind VAO
        glUseProgram(0);  // Unbind shader program
    }
}

