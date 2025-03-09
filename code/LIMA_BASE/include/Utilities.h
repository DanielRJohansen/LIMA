#pragma once

#include "LimaTypes.cuh"

#include <assert.h>
#include <cmath>
#include <concepts>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <ranges>
#include <string>
#include <vector>


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

    static void genericErrorCheckNoSync(const char* text) {
        if constexpr (SYNC_ALL_KERNELS)
            cudaDeviceSynchronize();

        cudaError_t cuda_status = cudaGetLastError();
        if (cuda_status != cudaSuccess) {
            std::cout << "\nCuda error code: " << cuda_status << " - " << cudaGetErrorString(cuda_status) << std::endl;
            fprintf(stderr, text);
            throw std::runtime_error("genericErrorCheck failed");
        }
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
    LimaLogger(const LogMode mode, EnvMode envmode, const std::string& name, const std::filesystem::path& workfolder=""); // With no workfolder, the logger simply wont putput anything to a file
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

// Lima Algorithm Library
namespace LAL {

    // TODO: Make a lowest-level file for these type agnostic algo's so we can use them in limatypes.cuh
    //template <typename T>
    //__device__ __host__ static T max(const T l, const T r) {
    //    return r > l ? r : l;
    //}
    
    bool constexpr isPowerOf2(int n) {
        return (n > 0) && ((n & (n - 1)) == 0);
    }


    constexpr int powi(int base, int exp) {
        int res = 1;
        for (int i = 0; i < exp; i++) {
			res *= base;
		}
        return res;
    }



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
