#pragma once

#include <iostream>
#include <unordered_map>
#include <vector>
#include <functional>
#include <algorithm>
#include <filesystem>
#include <format>

namespace fs = std::filesystem;

class ArgParser {
public:
    using Callback = std::function<void(const std::vector<std::string>& args)>;
    using FlagCallback = std::function<void()>;

    ArgParser(const std::string& helpText)
        : helpText(helpText)
    {
        aliasMap.insert({ "-h",     "-help" });
        aliasMap.insert({ "-Help",  "-help" });
    }

    void AddOption(const std::vector<std::string>& aliases, bool required, Callback callback) {
        for (auto& a : aliases) aliasMap[a] = aliases[0];
        options[aliases[0]] = { required, std::move(callback) };
    }

    void AddOption(const std::vector<std::string>& aliases, bool required, std::string& target) {
        AddOption(aliases, required,
            [&target, name = aliases[0]](auto const& args) {
                target = args[0];
            });
    }

    void AddOption(const std::vector<std::string>& aliases, bool required, fs::path& target) {
        AddOption(aliases, required,
            [&target](auto const& args) {
                target = args[0];
            });
    }

    void AddOption(const std::vector<std::string>& aliases, bool required, int& target) {
        AddOption(aliases, required,
            [&target, name = aliases[0]](auto const& args) {
                try {
                    target = std::stoi(args[0]);
                }
                catch (...) {
                    std::cerr << std::format("Argument {} expected an integer, got: {}\n", name, args[0]);
                    std::exit(1);
                }
            });
    }

    void AddOption(const std::vector<std::string>& aliases, bool required, float& target) {
        AddOption(aliases, required,
            [&target, name = aliases[0]](auto const& args) {
                try {
                    target = std::stof(args[0]);
                }
                catch (...) {
                    std::cerr << std::format("Argument {} expected a float, got: {}\n", name, args[0]);
                    std::exit(1);
                }
            });
    }
    void AddOption(const std::vector<std::string>& aliases, bool required, std::optional<float>& target) {
        AddOption(aliases, required,
            [&target, name = aliases[0]](auto const& args) {
                try {
                    target = std::stof(args[0]);
                }
                catch (...) {
                    std::cerr << std::format("Argument {} expected a float, got: {}\n", name, args[0]);
                    std::exit(1);
                }
            });
    }
    void AddOption(const std::vector<std::string>& aliases, bool required, Float3& target, bool allowInitFrom1=false) {
        AddOption(aliases, required,
            [&target, name = aliases[0], allowInitFrom1](auto const& args) {
                if (args.size() == 3 || (allowInitFrom1 && args.size() == 1)) {}
                else {
                    std::cerr << std::format("Argument {} expected 3 floats, got {}\n", name, args.size());
                    std::exit(1);
                }

                try {
                    if (args.size() == 3) {
                        target.x = std::stof(args[0]);
                        target.y = std::stof(args[1]);
                        target.z = std::stof(args[2]);
                    }
                    else
						target = Float3{ std::stof(args[0])};
                }
                catch (...) {
                    std::cerr << std::format("Argument {} expected 3 floats, invalid value in {}\n", name, args[0]);
                    std::exit(1);
                }
            });
    }

    void AddFlag(const std::vector<std::string>& aliases, FlagCallback callback) {
        for (auto& a : aliases) aliasMap[a] = aliases[0];
        flags[aliases[0]] = std::move(callback);
    }

    void Parse(int argc, char** argv) {
        std::unordered_map<std::string, bool> seen;
        for (int i = 2; i < argc; ++i) {
            std::string key = argv[i];
            if (!aliasMap.contains(key)) {
                std::cerr << std::format("Unknown argument: {}\n", key);
                std::cout << helpText;
                std::exit(1);
            }

            auto const& name = aliasMap[key];
            if (name == "-help") {
                std::cout << helpText;
                std::exit(0);
            }

            if (auto oit = options.find(name); oit != options.end()) {
                seen[name] = true;
                std::vector<std::string> args;
                int j = i + 1;
                while (j < argc && argv[j][0] != '-') {
                    args.emplace_back(argv[j++]);
                }
                if (args.empty()) {
                    std::cerr << std::format("Argument {} expected a value\n", name);
                    std::exit(1);
                }
                oit->second.callback(args);
                i = j - 1;
            }
            else if (auto fit = flags.find(name); fit != flags.end()) {
                fit->second();
            }
        }

        for (auto& [name, opt] : options) {
            if (opt.required && !seen[name]) {
                std::cerr << std::format("Missing required option: {}\n", name);
                std::cout << helpText;
                std::exit(1);
            }
        }
    }

private:
    struct Option { bool required; Callback callback; };
    std::unordered_map<std::string, std::string>      aliasMap;
    std::unordered_map<std::string, Option>           options;
    std::unordered_map<std::string, FlagCallback>     flags;
    std::string                                       helpText;
};
