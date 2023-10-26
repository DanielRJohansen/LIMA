#include "Utilities.h"
#include <filesystem>
#include <iostream>
#include "Printer.h"

using namespace LIMA_Print;
using std::string;

LimaLogger::LimaLogger(LogMode lm, EnvMode em, const string& name, const string& workfolder)
    : logmode(lm)
    , envmode{em}
    , enable_logging(workfolder !="")
    , log_dir(workfolder + "/logs/")
{
    if (enable_logging) {
        std::filesystem::create_directories(log_dir);
        const auto logFilePath = log_dir + name + "_log.txt";
        logFile.open(logFilePath, std::ios::out);
    }    
}

LimaLogger::~LimaLogger() {
    if (enable_logging && logFile.is_open()) {
        logFile.close();
    }
}

void LimaLogger::startSection(const std::string& input)
{
    if (envmode != Headless) {
        printH2(input, true, false);
    }

    logToFile(input);
}

void LimaLogger::print(const std::string& input, const bool log) 
{
    if (envmode != Headless) {
        string input_copy = input;
        if (logmode == compact) {
            if (clear_next) { clearLine(); }
            if (!input_copy.empty() && input_copy.back() == '\n') {
                input_copy.back() = ' ';
                clear_next = true;  // Only clear of this was end of line
            }
        }
        std::cout << input_copy;
    }


    if (log) {
        logToFile(input);
    }
}

void LimaLogger::finishSection(const string& str) {
    if (envmode != Headless) {
        if (logmode == compact) {
            std::cout << "\n";
        }

        printH2(str, false, true);
    }
    
    logToFile(str);
}

void LimaLogger::logToFile(const string& str)
{
	if (!enable_logging) { return; }
	if (!logFile.is_open()) {
		std::cerr << "Error: Log file is not open" << std::endl;
		return;
	}

	logFile << str << std::endl;
}

void LimaLogger::clearLine() {
    //int a = static_cast<int>(std::cerr.tellp()) - 1;
    //std::cout << "\r" << std::flush;
    //std::cout << std::string(std::max(0, static_cast<int>(std::cerr.tellp()) - 1), ' ') << "\r" << std::flush;
    std::cout << "\033[2K\r" << std::flush;
    clear_next = false;
}