#pragma once

#include <glm.hpp>
#include <GL/glew.h>
#include <gtc/matrix_transform.hpp>
#include <gtc/type_ptr.hpp>

#include <string>
#include <iostream>



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

    void Expand(size_t byteSize) {
        if (byteSize > currentSize) {
            Resize(byteSize);
        }
    }
    const size_t Capacity() const {
        return currentSize;
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

    template<typename T>
    std::vector<T> GetData() {
        std::vector<T> out(currentSize / sizeof(T));
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, bufferID);
        void* mappedData = glMapBuffer(GL_SHADER_STORAGE_BUFFER, GL_READ_ONLY);
        if (mappedData) {
            // Copy the data from the mapped buffer to the vector
            std::memcpy(out.data(), mappedData, currentSize);
            // Unmap the buffer
            glUnmapBuffer(GL_SHADER_STORAGE_BUFFER);
        }
        else {
            // Handle error (e.g., throw an exception or log a message)
            throw std::runtime_error("Failed to map OpenGL buffer for reading.");
        }

        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
        return out;
    }

    GLuint GetID() const {
        return bufferID;
    }
};
