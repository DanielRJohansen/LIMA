#pragma once

#include <assert.h>
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include "LimaTypes.cuh"
#include <filesystem>

#include <memory>
namespace fs = std::filesystem;

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
            throw std::runtime_error("genericErrorCheck failed");
        }
    }

    static void genericErrorCheckNosync(const char* text) {
        if constexpr (!LIMA_PUSH)
            cudaDeviceSynchronize();

        cudaError_t cuda_status = cudaGetLastError();
        if (cuda_status != cudaSuccess) {
            std::cout << "\nCuda error code: " << cuda_status << " - " << cudaGetErrorString(cuda_status) << std::endl;
            fprintf(stderr, text);
            throw std::runtime_error("genericErrorCheck failed");
        }
    }
    static void genericErrorCheckDebugOnly(const char* text) {
#ifdef _DEBUG
        cudaDeviceSynchronize();
        cudaError_t cuda_status = cudaGetLastError();
        if (cuda_status != cudaSuccess) {
            std::cout << "\nCuda error code: " << cuda_status << " - " << cudaGetErrorString(cuda_status) << std::endl;
            fprintf(stderr, text);
            throw std::runtime_error("genericErrorCheck failed");
        }
#endif
    }

    static void genericErrorCheck(const cudaError_t cuda_status) {
        if (cuda_status != cudaSuccess) {
            std::cout << "\nCuda error code: " << cuda_status << " - " << cudaGetErrorString(cuda_status) << std::endl;
            throw std::runtime_error("genericErrorCheck failed");
        }
    }
}




class LimaLogger {
public:
    enum LogMode {
        normal,
        compact
    };   
    LimaLogger() {}
    LimaLogger(const LimaLogger&) = delete;
    LimaLogger(const LogMode mode, EnvMode envmode, const std::string& name, const fs::path& workfolder=""); // With no workfolder, the logger simply wont putput anything to a file
    ~LimaLogger();

    void startSection(const std::string& input);
    void print(const std::string& input, bool log=true);
    void finishSection(const std::string& str);
    
    template <typename T>
    void printToFile(const std::string& filename, const std::vector<T>& data) const {
        // Does nothing
    }


private:
    LogMode logmode{};
    EnvMode envmode{};
    //std::string logFilePath;
    const std::string log_dir;
    std::ofstream logFile;
    const bool enable_logging{ false };

    void logToFile(const std::string& str);
    void clearLine();
    bool clear_next = false;
};

static std::unique_ptr<LimaLogger> makeLimaloggerBareboned(const std::string& name) {
    return std::make_unique<LimaLogger>(LimaLogger::LogMode::compact, EnvMode::Headless, name);
}
#include <cmath>
#include <concepts>
#include <ranges>
// Lima Algorithm Library
namespace LAL {

    // TODO: Make a lowest-level file for these type agnostic algo's so we can use them in limatypes.cuh
    //template <typename T>
    //__device__ __host__ static T max(const T l, const T r) {
    //    return r > l ? r : l;
    //}
    


    constexpr int powi(int base, int exp) {
        int res = 1;
        for (int i = 0; i < exp; i++) {
			res *= base;
		}
        return res;
    }
    constexpr int ceilFloatToInt(float num) {
        return (static_cast<float>(static_cast<int>(num)) == num)
            ? static_cast<int>(num)
            : static_cast<int>(num) + ((num > 0) ? 1 : 0);
    }

    constexpr float pow_constexpr(float base, int exp) {
        return (exp == 0) ? 1 : base * pow_constexpr(base, exp - 1);
    }

    // This is not sustainable....
    constexpr float log2f(float val) {
        return (val < 2.0f) 
            ? (val - 1.0f) / (val + 1.0f) + (1.0f / 3.0f) 
                * pow_constexpr((val - 1.0f) / (val + 1.0f), 3) + (1.0f / 5.0f) 
                * pow_constexpr((val - 1.0f) / (val + 1.0f), 5) 
            : 1.0f + log2f(val / 2.0f);
    }

    template<typename T>
    class optional {
    private:
        bool _hasValue = false;
        T _value;

    public:
        optional() : _hasValue(false), _value() {}

        optional(const T& value) : _hasValue(true), _value(value) {}

        optional(T&& value) : _hasValue(true), _value(std::move(value)) {}

        optional& operator=(const T& value) {
            this->_value = value;
            _hasValue = true;
            return *this;
        }

        optional& operator=(T&& value) {
            this->_value = std::move(value);
            _hasValue = true;
            return *this;
        }

        T value() const {
            if (!_hasValue) {
                throw std::runtime_error("No value present in optional.");
            }
            return _value;
        }

        operator bool() const {
            return _hasValue;
        }

        bool hasValue() const { return _hasValue; }
    };    


    //float LargestDiff(const Float3 queryPoint, const std::span<Float3>& points);

    template <std::ranges::range ContainerType>
    float LargestDiff(const Float3 queryPoint, const ContainerType& points) {
        float maxDiff = 0.f;
        for (const Float3& p : points) {
            maxDiff = std::max(maxDiff, (p - queryPoint).len()); // Compute lenSq if optimizing
        }
        return maxDiff;
    }
}
