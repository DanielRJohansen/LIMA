#pragma once

#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <unistd.h>

// Are these necessary?
#include <sys/stat.h>
#include <sys/types.h>

namespace fs = std::filesystem;

namespace SelfRecompile {
    struct UserConstantInfo {
        std::string type;
        std::string value;
    };


    std::map<std::string, UserConstantInfo> readDefaultConstants(const std::string& filename) {
        std::ifstream infile(filename);
        std::string line;
        std::map<std::string, UserConstantInfo> defaultConstants;

        while (std::getline(infile, line))
        {
            // Check if the line starts with a "/" (comment)
            if (!line.empty() && line[0] == '/') {
                continue; // Ignore comments
            }

            std::istringstream iss(line);
            std::string type1, type2, key, equals, value;

            if (iss >> type1 >> type2 >> key >> equals >> value) {
                UserConstantInfo info = { type1 + " " + type2, value };
                // Remove trailing semicolon if it exists
                if (value.back() == ';') {
                    value.pop_back();
                }
                defaultConstants[key] = info;
            }
        }
        return defaultConstants;
    }

    void readAndOverrideConstants(const std::string& filename, std::map<std::string, UserConstantInfo>& constants) {
        std::ifstream infile(filename);
        std::string line;

        while (std::getline(infile, line)) {
            std::istringstream iss(line);
            std::string key, value;
            if (std::getline(iss, key, '=') && std::getline(iss, value)) {
                //std::cout << std::format("Key {} value {}\n", key, value);

                if (constants.find(key) != constants.end()) {
                    constants[key].value = value;
                }
            }
        }
    }

    void writeConstantsToFile(const std::string& filename, const std::map<std::string, UserConstantInfo>& constants) {
        std::ofstream outfile(filename);
        outfile << "#pragma once\n\n";
        for (const auto& pair : constants) {
            outfile << pair.second.type << " " << pair.first << " = " << pair.second.value << ";\n";
        }
    }

    void overrideUserParams() {
        const std::string currentDirectory = fs::current_path().string();
        std::map<std::string, UserConstantInfo> constants = readDefaultConstants("/opt/LIMA/code/LIMA_BASE/include/DefaultUserConstants.h");

        const std::string params_path = currentDirectory + "/sim_params.txt";
        if (!fs::exists(params_path)) {
            throw std::runtime_error(std::format("Expected file {}, but it was not found", params_path).c_str());
        }
        readAndOverrideConstants(params_path, constants);

        writeConstantsToFile("UserConstants.h", constants);
        if (!fs::exists("UserConstants.h")) {
            throw std::runtime_error(std::format("Expected file {}, but it was not found", params_path).c_str());
        }
    }




















    void clearDirectory(const std::string& path, bool keep_directory) {
        if (fs::exists(path) && fs::is_directory(path)) {
            std::cout << "clearing dir" << path << "...";
            fs::remove_all(path);
            if (keep_directory)
                fs::create_directory(path);
            std::cout << "done!\n";
        }
    }

    void copyFiles(const std::string& src, const std::string& dest) {
        try {
            fs::copy(src, dest, fs::copy_options::overwrite_existing | fs::copy_options::recursive);
            std::cout << "Files copied from " << src << " to " << dest << std::endl;
        }
        catch (const fs::filesystem_error& e) {
            std::cerr << e.what() << std::endl;
        }
    }

    // Copies everything from /opt/LIMA to ~/LIMA
    void copySourceToUserProgram() {
        const std::string home = getenv("HOME");
        const std::string userprogramDir = home + "/LIMA";
        const std::string optDir = "/opt/LIMA";

        fs::create_directories(userprogramDir);

        copyFiles(optDir, userprogramDir);
    }

    bool runSystemCommand(const std::string& command, const std::string& logFile) {
        std::string fullCommand = command + " > " + logFile + " 2>&1";
        int result = std::system(fullCommand.c_str());
        if (result != 0) {
            std::cerr << "Command failed: " << command << std::endl;
            std::ifstream log(logFile);
            if (log.is_open()) {
                std::cerr << log.rdbuf();
                log.close();
            }
            return false;
        }
        return true;
    }

    int recompile() {
        const std::string home = getenv("HOME");
        const std::string programDir = home + "/LIMA/";
        const std::string buildDir = home + "/LIMA/build";
        //const std::string sourceDir = home + "/LIMA/source";
        const std::string applicationsDir = home + "/LIMA/applications";
        const std::string logFile = buildDir + "/limabuild.log";

        // Copy UserConstants.h
        fs::copy("UserConstants.h", programDir + "code/LIMA_BASE/include/UserConstants.h", fs::copy_options::overwrite_existing);

        fs::create_directories(buildDir);
        fs::create_directories(applicationsDir);
        clearDirectory(buildDir, true);

         // Change current_path to buildDir, then compile files, and finally switch back to the working dir
        fs::path work_dir = fs::current_path();
        fs::current_path(buildDir);
        // Run cmake and make
        if (!runSystemCommand("cmake " + programDir + " -Wno-dev", logFile) ||
            !runSystemCommand("make install -j", logFile)) {
            fs::current_path(work_dir);
            return 1;
        }
        // Move LIMA_TESTS/limatests
        fs::rename("code/LIMA_TESTS/limatests", "../applications/limatests");
        fs::current_path(work_dir); // Change back to previous path
        clearDirectory(buildDir, false);
        
        return 0;
    }




    int autoRecompile() 
    {
        copySourceToUserProgram();

        // Override default userparams with sim_params from user
        overrideUserParams();

        std::printf("Optimizing LIMA engine for your simulation parameters (This should take approx 1 minute)\n");
        return recompile();
    }
}
