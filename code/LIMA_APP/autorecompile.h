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
        char cwd[1024];
        getcwd(cwd, sizeof(cwd));
        std::string currentDirectory(cwd);
        std::map<std::string, UserConstantInfo> constants = readDefaultConstants("/opt/LIMA/code/LIMA_BASE/include/DefaultUserConstants.h");

        const std::string params_path = currentDirectory + "/sim_params.txt";
        if (!fs::exists(params_path)) {
            throw std::runtime_error(std::format("Expected file {}, but it was not found", params_path).c_str());
        }
        readAndOverrideConstants(params_path, constants);

        writeConstantsToFile("./UserConstants.h", constants);
    }




















    void clearDirectory(const std::string& path) {
        if (fs::exists(path) && fs::is_directory(path)) {
            fs::remove_all(path);
            fs::create_directory(path);
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

    int copySourceToUserProgram() {
        const std::string home = getenv("HOME");
        const std::string userprogramDir = home + "/LIMA";
        const std::string optDir = "/opt/LIMA";

        fs::create_directories(userprogramDir);

        //clearDirectory(userprogramDir + "/source/code");  // In case there was already code there

        copyFiles(optDir, userprogramDir);

        return 0;
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
        const std::string applicationsDir = home + "/LIMA/applications";
        const std::string logFile = buildDir + "/limabuild.log";

        // Copy UserConstants.h
        fs::copy("UserConstants.h", programDir + "code/LIMA_BASE/include/UserConstants.h", fs::copy_options::overwrite_existing);

        fs::create_directories(buildDir);
        fs::create_directories(applicationsDir);
        clearDirectory(buildDir);

        // Change to build directory
        fs::current_path(buildDir);

        // Run cmake and make
        if (!runSystemCommand("cmake " + programDir + " -Wno-dev", logFile) ||
            !runSystemCommand("make install", logFile)) {
            return 1;
        }

        // Move LIMA_TESTS/limatests
        fs::rename("LIMA_TESTS/limatests", "../limatests");

        return 0;
    }



    int autoRecompile() 
    {
        //const int err = system((source_dir + "copyToUserProgram.sh").c_str());
        //if (err) return err;
        // Copy all code to to ~/LIMA/source
        copySourceToUserProgram();

        // Override default userparams with sim_params from user
        overrideUserParams();

        // Call compile script
        std::printf("Optimization LIMA engine for your simulation parameters (This should take approx 1 minute)\n");
        // This call goes to the wrong dir, but the script will cd to the right folder before it compiles
        //return system((source_dir + "recompile.sh").c_str());
        recompile();
    }
}
