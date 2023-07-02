#pragma once

#include <assert.h>
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include "LimaTypes.cuh"

#include <device_launch_parameters.h>
#include <cuda_runtime_api.h>

namespace LIMA_UTILS {

	static int roundUp(int numToRound, int multiple)
	{
		assert(multiple);
		return ((numToRound + multiple - 1) / multiple) * multiple;
	}

    static void genericErrorCheck(const char* text) {
        cudaDeviceSynchronize();
        cudaError_t cuda_status = cudaGetLastError();
        if (cuda_status != cudaSuccess) {
            std::cout << "\nCuda error code: " << cuda_status << " - " << cudaGetErrorString(cuda_status) << std::endl;
            fprintf(stderr, text);
            exit(1);
        }
    }
    static void genericErrorCheck(const cudaError_t cuda_status) {
        if (cuda_status != cudaSuccess) {
            std::cout << "\nCuda error code: " << cuda_status << " - " << cudaGetErrorString(cuda_status) << std::endl;
            exit(1);
        }
    }
}




class LimaLogger {
public:
    enum LogMode {
        normal,
        compact
    };   

    LimaLogger(const LimaLogger&) = delete;
    LimaLogger(const LogMode mode, EnvMode envmode, const std::string& name, const std::string& workfolder="");
    ~LimaLogger();

    void startSection(const std::string& input);
    void print(const std::string& input, bool log=true);
    void finishSection(const std::string& str);
    
    template <typename T>
    void printToFile(const std::string& filename, const std::vector<T>& data) const {
        FILE* file;

        const std::string filepath = (log_dir + filename);
        if (!fopen_s(&file, filepath.c_str(), "wb")) {
            fwrite(data.data(), sizeof(T), data.size(), file);
            fclose(file);
        }
    }


private:
    LogMode logmode;
    EnvMode envmode;
    //std::string logFilePath;
    const std::string log_dir;
    std::ofstream logFile;
    const bool enable_logging{ false };

    void logToFile(const std::string& str);
    void clearLine();
    bool clear_next = false;
};

// Lima Algorithm Library
namespace LAL {


    template<typename T>
    class optional {
    private:
        bool hasValue = false;
        T _value;

    public:
        optional() : hasValue(false), _value() {}

        optional(const T& value) : hasValue(true), _value(value) {}

        optional(T&& value) : hasValue(true), _value(std::move(value)) {}

        optional& operator=(const T& value) {
            this->_value = value;
            hasValue = true;
            return *this;
        }

        optional& operator=(T&& value) {
            this->_value = std::move(value);
            hasValue = true;
            return *this;
        }

        T value() const {
            if (!hasValue) {
                throw std::exception("No value present in optional.");
            }
            return _value;
        }

        operator bool() const {
            return hasValue;
        }
    };    
}