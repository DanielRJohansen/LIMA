#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <string>
#include <unistd.h>
#include <stdlib.h>

namespace SelfRecompile {
    struct UserConstantInfo {
        std::string type;
        std::string value;
    };


    std::map<std::string, ConstantInfo> readDefaultConstants(const std::string& filename) {
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

    void readAndOverrideConstants(const std::string& filename, std::map<std::string, ConstantInfo>& constants) {
        std::ifstream infile(filename);
        std::string line;

        while (std::getline(infile, line)) {
            std::istringstream iss(line);
            std::string key, value;
            if (std::getline(iss, key, '=') && std::getline(iss, value)) {
                if (constants.find(key) != constants.end()) {
                    constants[key].value = value;
                }
            }
        }
    }

    void writeConstantsToFile(const std::string& filename, const std::map<std::string, ConstantInfo>& constants) {
        std::ofstream outfile(filename);

        for (const auto& pair : constants) {
            outfile << pair.second.type << " " << pair.first << " = " << pair.second.value << ";\n";
        }
    }

    void overrideUserParams() {
        char cwd[1024];
        getcwd(cwd, sizeof(cwd));
        std::string currentDirectory(cwd);

        std::map<std::string, ConstantInfo> constants = readDefaultConstants("/home/opt/LIMA/Applications/LIMA_BASE/src/DefaultUserConstants.h");

        readAndOverrideConstants(currentDirectory + "/simparams.txt", constants);

        writeConstantsToFile("/home/opt/LIMA/Applications/LIMA_BASE/src/UserConstants.h", constants);
    }

    int autoRecompile() 
    {
        overrideUserParams();

        // Call compile script
        system("./home/opt/LIMA/Applications/recompile.sh");
    }
}
