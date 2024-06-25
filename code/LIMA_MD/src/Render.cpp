#include "DisplayV2.h"

#include <glm.hpp>
#include <gtc/matrix_transform.hpp>
#include <gtc/type_ptr.hpp>





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
    if (!drawBoxShaderProgram.has_value()) {
        drawBoxShaderProgram.emplace(0);
        SetupDrawBoxPipeline(*drawBoxShaderProgram);
    }

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








GLuint drawAtomsVBO, drawAtomsSSBO;

void SetupDrawAtomsPipeline(size_t numAtoms, GLuint& drawAtomsShaderProgram) {
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
    glGenBuffers(1, &drawAtomsSSBO);

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




void Display::DrawAtoms(const std::vector<RenderAtom>& atoms) {
    const int numVerticesPerAtom = 12; // Number of vertices to define a disc (30 for the circle + 1 center + 1 to close the loop)

    if (!drawAtomsShaderProgram.has_value()) {
        drawAtomsShaderProgram.emplace(0);
        SetupDrawAtomsPipeline(atoms.size() * numVerticesPerAtom, *drawAtomsShaderProgram);
    }

    glUseProgram(*drawAtomsShaderProgram);

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, drawAtomsSSBO);
    glBufferData(GL_SHADER_STORAGE_BUFFER, atoms.size() * sizeof(RenderAtom), atoms.data(), GL_DYNAMIC_DRAW);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, drawAtomsSSBO);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

    const GLuint mvpLocation = glGetUniformLocation(*drawAtomsShaderProgram, "MVP");
    const glm::mat4 MVP = GetMVPMatrix(camera_distance, camera_pitch * rad2deg, camera_yaw * rad2deg, screenWidth, screenHeight);
    glUniformMatrix4fv(mvpLocation, 1, GL_FALSE, glm::value_ptr(MVP));

    const GLuint numAtomsLocation = glGetUniformLocation(*drawAtomsShaderProgram, "numAtoms");
    glUniform1i(numAtomsLocation, atoms.size());

    const GLuint numVerticesPerAtomLocation = glGetUniformLocation(*drawAtomsShaderProgram, "numVerticesPerAtom");
    glUniform1i(numVerticesPerAtomLocation, numVerticesPerAtom);

    glDrawArraysInstanced(GL_TRIANGLE_FAN, 0, numVerticesPerAtom, atoms.size());

    glBindVertexArray(0);
    glUseProgram(0);
}


























void Display::TerminateGLEW() {
    // Cleanup
    glDeleteVertexArrays(1, &drawAtomsVBO);
    glDeleteVertexArrays(1, &boxVAO);
    glDeleteVertexArrays(1, &boxVBO);
    glDeleteVertexArrays(1, &drawAtomsSSBO);
    

    if (drawAtomsShaderProgram.has_value())
        glDeleteProgram(*drawAtomsShaderProgram);
    if (drawBoxShaderProgram.has_value())
        glDeleteProgram(*drawBoxShaderProgram);
}