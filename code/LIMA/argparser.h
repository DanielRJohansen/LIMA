#pragma once

#include <iostream>
#include <unordered_map>
#include <vector>
#include <functional>
#include <algorithm>
#include <filesystem>

namespace fs = std::filesystem;

class ArgParser {
public:
    using Callback = std::function<void(int& indexOfOption, char** argv, int argc)>;
    using FlagCallback = std::function<void()>;

    ArgParser(const std::string& helpText) : helpText(helpText) {
        aliasMap.insert({"-h", "-help"});
        aliasMap.insert({"-Help", "-help"});
    }

    void AddOption(const std::vector<std::string>& aliases, bool required, Callback callback) {
        for (const auto& alias : aliases) {
            aliasMap.insert({alias, aliases[0]});
        }
        options[aliases[0]] = { required, callback };
    }

    void AddOption(const std::vector<std::string>& aliases, bool required, std::string& target) {
        AddOption(aliases, required, [&target](int indexOfOption, char** argv, int) {
            target = argv[indexOfOption];
            });
    }

    void AddOption(const std::vector<std::string>& aliases, bool required, std::filesystem::path& target) {
        AddOption(aliases, required, [&target](int indexOfOption, char** argv, int) {
            target = argv[indexOfOption];
            });
    }

    void AddOption(const std::vector<std::string>& aliases, bool required, int& target) {
        const std::string argMainName = aliases[0];
        AddOption(aliases, required, [&target, argMainName](int indexOfOption, char** argv, int) {
            try {
                target = std::stoi(argv[indexOfOption]);
            }
            catch (std::invalid_argument) {
                std::cerr << std::string("Argument ") + argMainName + "expected an integer, got value: " + argv[indexOfOption] + "\n";
                exit(1);
            }
            });
    }

    void AddOption(const std::vector<std::string>& aliases, bool required, Float3& target) {
        const std::string argMainName = aliases[0];
        AddOption(aliases, required, [&target, argMainName](int indexOfOption, char** argv, int argc) {
            if (indexOfOption+2 >= argc) {
				std::cerr << std::string("Argument ") + argMainName + "expected 3 floats, got " + std::to_string(argc-indexOfOption+1) + "\n";
				exit(1);
			}

            try {
                target.x = std::stof(argv[indexOfOption]);
                target.y = std::stof(argv[++indexOfOption]);
                target.z = std::stof(argv[++indexOfOption]);
            }
            catch (std::invalid_argument) {
                std::cerr << std::string("Argument ") + argMainName + "expected an 3 floats, got value: " + argv[indexOfOption] + "\n";
                exit(1);
            }
            });
    }

    void AddFlag(const std::vector<std::string>& aliases, FlagCallback callback) {
        for (const auto& alias : aliases) {
            aliasMap.insert({alias, aliases[0]});
        }
        flags[aliases[0]] = callback;
    }

    void Parse(int argc, char** argv) {
        std::unordered_map<std::string, bool> matchedOptions;
        for (int i = 2; i < argc; ++i) {
            const std::string argAlias = argv[i];
            if (!aliasMap.contains(argAlias)) {
                std::cerr << "Unknown argument: " << argAlias << "\n";
                std::cout << helpText;
                exit(1);
            }
            const std::string arg = aliasMap[argAlias];

            if (arg == "-help") {
                std::cout << helpText;
                exit(0);
            }


            auto optIt = options.find(arg);
            auto flagIt = flags.find(arg);
            if (optIt != options.end()) {
                matchedOptions[arg] = true;
                if (i + 1 < argc && argv[i + 1][0] != '-') {
                    //optIt->second.callback(argv[++i]);
                    optIt->second.callback(++i, argv, argc);
                }
                else {
                    std::cerr << "Argument " + arg + " expected a value\n";
                    exit(1);
                }
            }
            else if (flagIt != flags.end()) {
                flagIt->second();
            }
            else {
                std::cerr << "Unknown error encountered";
            }
        }

        // Check for required options
        for (const auto& [alias, option] : options) {
            if (option.required && !matchedOptions[alias]) {
                std::cerr << "Missing required option: " << alias << "\n";
                std::cout << helpText;
                exit(1);
            }
        }
    }
      
private:
    struct Option {
        bool required;
        Callback callback;
    };

    std::unordered_map<std::string, std::string> aliasMap;

    std::unordered_map<std::string, Option> options;
    std::unordered_map<std::string, FlagCallback> flags;
    std::string helpText;
};
