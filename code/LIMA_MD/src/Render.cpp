#include "DisplayV2.h"

#include <glm.hpp>
#include <gtc/matrix_transform.hpp>
#include <gtc/type_ptr.hpp>


#include <cuda_gl_interop.h>

#include <GLShader.h>

// ------------------------------------------------------------- Render functions used during sim ------------------------------------------------------------- //

static glm::mat4 GetMVPMatrix(float camera_distance, float camera_pitch, float camera_yaw, int screenWidth, int screenHeight) {
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
























void Display::initializePipeline(size_t numAtoms) {
    drawAtomsShaderProgram.emplace(0);

	// Set up the draw atoms pipeline
	SetupDrawAtomsPipeline(numAtoms, *drawAtomsShaderProgram, renderAtomsBuffer, &renderAtomsBufferCudaResource);
    

    glEnable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}






void Display::TerminateGLEW() {
    if (!drawAtomsShaderProgram.has_value())
        return; // Pipeline was never initialized, so nothing to clean up

    // Cleanup
    glDeleteVertexArrays(1, &drawAtomsVBO);
    glDeleteVertexArrays(1, &renderAtomsBuffer.value());
    

    if (drawAtomsShaderProgram.has_value())
        glDeleteProgram(*drawAtomsShaderProgram);
    if (drawBoxShaderProgram.has_value())
        glDeleteProgram(*drawBoxShaderProgram);
}







