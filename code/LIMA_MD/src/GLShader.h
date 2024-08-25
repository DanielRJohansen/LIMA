#pragma once

#include "Utilities.h"

#include <glm.hpp>
#include <GL/glew.h>
#include <gtc/matrix_transform.hpp>
#include <gtc/type_ptr.hpp>

#include <string>
#include <iostream>


class Shader {
public:
    GLuint programID;

    Shader(const std::string& vertexCode, const std::string& fragmentCode) {
        GLuint vertexShader = compileShader(vertexCode, GL_VERTEX_SHADER);
        GLuint fragmentShader = compileShader(fragmentCode, GL_FRAGMENT_SHADER);
        programID = linkProgram(vertexShader, fragmentShader);
    }
    ~Shader() {
		glDeleteProgram(programID);
	}

    void use() {
        glUseProgram(programID);
    }

    void SetUniformI(const std::string& name, int value) {
        glUniform1i(glGetUniformLocation(programID, name.c_str()), value);
    }

    //void setUniform(const std::string& name, float value) {
    //    glUniform1f(glGetUniformLocation(programID, name.c_str()), value);
    //}

    //void setUniform(const std::string& name, const glm::vec2& value) {
    //    glUniform2f(glGetUniformLocation(programID, name.c_str()), value.x, value.y);
    //}

    void SetUniformFloat3(const std::string& name, const Float3& value) {
        glUniform3f(glGetUniformLocation(programID, name.c_str()), value.x, value.y, value.z);
    }

    void SetUniformFloat4(const std::string& name, const float4& value) {
        glUniform4f(glGetUniformLocation(programID, name.c_str()), value.x, value.y, value.z, value.w);
    }

    void SetUniformMat4(const std::string& name, const glm::mat4& mat) {
        glUniformMatrix4fv(glGetUniformLocation(programID, name.c_str()), 1, GL_FALSE, glm::value_ptr(mat));
    }

private:
    GLuint compileShader(const std::string& source, GLenum shaderType) {
        GLuint shader = glCreateShader(shaderType);
        const char* src = source.c_str();
        glShaderSource(shader, 1, &src, nullptr);
        glCompileShader(shader);

        // Check for compilation errors
        GLint success;
        glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
        if (!success) {
            char infoLog[512];
            glGetShaderInfoLog(shader, 512, nullptr, infoLog);
            std::cerr << "ERROR::SHADER::COMPILATION_FAILED\n" << infoLog << std::endl;
            glDeleteShader(shader);
            throw std::runtime_error("Shader compilation failed");
        }

        return shader;
    }

    GLuint linkProgram(GLuint vertexShader, GLuint fragmentShader) {
        GLuint program = glCreateProgram();
        glAttachShader(program, vertexShader);
        glAttachShader(program, fragmentShader);
        glLinkProgram(program);

        // Check for linking errors
        GLint success;
        glGetProgramiv(program, GL_LINK_STATUS, &success);
        if (!success) {
            char infoLog[512];
            glGetProgramInfoLog(program, 512, nullptr, infoLog);
            std::cerr << "ERROR::PROGRAM::LINKING_FAILED\n" << infoLog << std::endl;
            glDeleteProgram(program);
            throw std::runtime_error("Program linking failed");
        }

        // Clean up shaders (they are no longer needed after linking)
        glDeleteShader(vertexShader);
        glDeleteShader(fragmentShader);

        return program;
    }
};

class SSBO {
private:
    GLuint bufferID;
    size_t currentSize;

public:
    SSBO() : currentSize(0) {
        glGenBuffers(1, &bufferID);
    }

    ~SSBO() {
        glDeleteBuffers(1, &bufferID);
    }

    void Bind(GLuint bindingIndex) const {
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, bindingIndex, bufferID);
    }

    void Unbind() const {
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
    }

    void Resize(size_t byteSize) {
        if (byteSize == currentSize) {
			return;
		}
		glBindBuffer(GL_SHADER_STORAGE_BUFFER, bufferID);
		glBufferData(GL_SHADER_STORAGE_BUFFER, byteSize, nullptr, GL_DYNAMIC_DRAW);
		currentSize = byteSize;
		glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
	}

    void SetData_FromCuda(const void* data, size_t byteSize) {
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, bufferID);
        if (byteSize > currentSize) {
            // Allocate new memory if the required size exceeds the current size.
            // glBufferData will automatically free the old memory if reallocating.
            glBufferData(GL_SHADER_STORAGE_BUFFER, byteSize, nullptr, GL_DYNAMIC_DRAW);
            currentSize = byteSize;
        }
        // Map buffer and copy data from CPU/GPU memory
        void* ptr = glMapBuffer(GL_SHADER_STORAGE_BUFFER, GL_WRITE_ONLY);
        cudaMemcpy(ptr, data, byteSize, cudaMemcpyDeviceToHost);
        glUnmapBuffer(GL_SHADER_STORAGE_BUFFER);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

        LIMA_UTILS::genericErrorCheck("SetData_FromCuda");
    }

    // Method to set data from a std::vector on the host
    template <typename T>
    void SetData(const std::vector<T>& dataVector) {
        size_t dataSize = dataVector.size() * sizeof(T);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, bufferID);
        if (dataSize > currentSize) {
            // Allocate new memory if the required size exceeds the current size.
            glBufferData(GL_SHADER_STORAGE_BUFFER, dataSize, nullptr, GL_DYNAMIC_DRAW);
            currentSize = dataSize;
        }
        // Map buffer and copy data from the host vector
        void* ptr = glMapBuffer(GL_SHADER_STORAGE_BUFFER, GL_WRITE_ONLY);
        memcpy(ptr, dataVector.data(), dataSize);
        glUnmapBuffer(GL_SHADER_STORAGE_BUFFER);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
    }

    GLuint GetID() const {
        return bufferID;
    }
};
