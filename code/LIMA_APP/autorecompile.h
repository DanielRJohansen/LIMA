#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <string>
#include <unistd.h>
#include <stdlib.h>
#include <filesystem>
#include <format>

namespace fs = std::filesystem;

namespace SelfRecompile {

    const std::string source_dir = "/opt/LIMA/source/";

    struct UserConstantInfo {
        std::string type;
        std::string value;
    };


    std::map<std::string, UserConstantInfo> readDefaultConstants(const std::string& filename) {
        std::ifstream infile(filename);
        std::string line;
        std::map<std::string, UserConstantInfo> defaultConstants;

        while (std::getline(infile, line)) {
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

        for (const auto& pair : constants) {
            outfile << pair.second.type << " " << pair.first << " = " << pair.second.value << ";\n";
        }
    }

    void overrideUserParams() {
        char cwd[1024];
        getcwd(cwd, sizeof(cwd));
        std::string currentDirectory(cwd);

        std::map<std::string, UserConstantInfo> constants = readDefaultConstants(source_dir + "LIMA_BASE/include/DefaultUserConstants.h");

        const std::string params_path = currentDirectory + "/sim_params.txt";
        if (!fs::exists(params_path)) {
            throw std::runtime_error(std::format("Expected file {}, but it was not found", params_path).c_str());
        }
        readAndOverrideConstants(params_path, constants);

        writeConstantsToFile("./UserConstants.h", constants);
    }

    int autoRecompile() 
    {
        const int err = system((source_dir + "copyToUserProgram.sh").c_str());
        if (err) return err;


        overrideUserParams();

        // Call compile script
        std::printf("Optimization LIMA engine for your simulation parameters (This should take approx 1 minute)\n");

        // This call goes to the wrong dir, but the script will cd to the right folder before it compiles
        return system((source_dir + "recompile.sh").c_str());
    }
}
