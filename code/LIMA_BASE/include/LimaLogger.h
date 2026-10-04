#pragma once

#include "LimaTypes.cuh"

#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

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
