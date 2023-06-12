#include "Utilities.h"
#include <filesystem>
#include <iostream>

LimaLogger::LimaLogger(const LogMode mode, const std::string& name, const std::string& workfolder, bool headless)
    : mode(mode)
    , enable_logging(workfolder !="")
    , log_dir(workfolder + "/logs/")
    , headless{headless}
{
    if (enable_logging) {
        std::filesystem::create_directories(log_dir);
        const auto logFilePath = log_dir + name + "_log.txt";
        logFile.open(logFilePath, std::ios::out | std::ios::app);
    }    
}

LimaLogger::~LimaLogger() {
    if (enable_logging && logFile.is_open()) {
        logFile.close();
    }
}

void LimaLogger::print(const std::string& input, const bool log) {
    // First log to file
    if (enable_logging) {
        if (!logFile.is_open()) {
            std::cerr << "Error: Log file is not open" << std::endl;
            return;
        }
        if (log) { logFile << input << std::endl; }
    }
    
    if (headless) { return; }

    // Then print to console
    auto input_copy = input;
    if (mode == compact) {
        if (clear_next) { clearLine(); }
        if (!input_copy.empty() && input_copy.back() == '\n') {
            input_copy.back() = ' ';
            clear_next = true;  // Only clear of this was end of line
        }
    }
    std::cout << input_copy;
}

void LimaLogger::finishSection() {
    if (mode == compact) {
        std::cout << "\n";
    }
}

void LimaLogger::clearLine() {
    //int a = static_cast<int>(std::cerr.tellp()) - 1;
    //std::cout << "\r" << std::flush;
    //std::cout << std::string(std::max(0, static_cast<int>(std::cerr.tellp()) - 1), ' ') << "\r" << std::flush;
    std::cout << "\033[2K\r" << std::flush;
    clear_next = false;
}