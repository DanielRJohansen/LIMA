#pragma once

#include <filesystem>
#include <format>
#include <functional>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace fs = std::filesystem;

class CliError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class HelpRequested {
public:
    explicit HelpRequested(std::string helpText) : helpText(std::move(helpText)) {}
    std::string helpText;
};

class ArgParser {
public:
    using Callback = std::function<void(const std::vector<std::string>&)>;
    using FlagCallback = std::function<void()>;

    explicit ArgParser(std::string helpText) : helpText(std::move(helpText)) {
        aliasMap.insert({ "--help", "--help" });
        aliasMap.insert({ "-h", "--help" });
        aliasMap.insert({ "-help", "--help" }); // Legacy spelling.
    }

    void AddOption(const std::vector<std::string>& aliases, bool required, Callback callback,
        std::size_t minValues = 1, std::size_t maxValues = (std::numeric_limits<std::size_t>::max)()) {
        ValidateAliases(aliases);
        for (const auto& alias : aliases) aliasMap[alias] = aliases[0];
        options[aliases[0]] = { required, minValues, maxValues, std::move(callback) };
    }

    void AddOption(const std::vector<std::string>& aliases, bool required, std::string& target) {
        AddOption(aliases, required, [&target](const auto& args) { target = args[0]; }, 1, 1);
    }

    void AddOption(const std::vector<std::string>& aliases, bool required, fs::path& target) {
        AddOption(aliases, required, [&target](const auto& args) { target = args[0]; }, 1, 1);
    }

    void AddOption(const std::vector<std::string>& aliases, bool required, int& target) {
        AddOption(aliases, required, [&target, name = aliases[0]](const auto& args) {
            try { target = std::stoi(args[0]); }
            catch (...) { throw CliError(std::format("option '{}' expects an integer; got '{}'", name, args[0])); }
        }, 1, 1);
    }

    void AddOption(const std::vector<std::string>& aliases, bool required, float& target) {
        AddOption(aliases, required, [&target, name = aliases[0]](const auto& args) {
            try { target = std::stof(args[0]); }
            catch (...) { throw CliError(std::format("option '{}' expects a number; got '{}'", name, args[0])); }
        }, 1, 1);
    }

    void AddOption(const std::vector<std::string>& aliases, bool required, std::optional<float>& target) {
        AddOption(aliases, required, [&target, name = aliases[0]](const auto& args) {
            try { target = std::stof(args[0]); }
            catch (...) { throw CliError(std::format("option '{}' expects a number; got '{}'", name, args[0])); }
        }, 1, 1);
    }

    void AddOption(const std::vector<std::string>& aliases, bool required, Float3& target, bool allowInitFrom1 = false) {
        AddOption(aliases, required, [&target, name = aliases[0], allowInitFrom1](const auto& args) {
            if (args.size() != 3 && !(allowInitFrom1 && args.size() == 1)) {
                throw CliError(std::format("option '{}' expects {} number(s); got {}",
                    name, allowInitFrom1 ? "one or three" : "three", args.size()));
            }
            try {
                if (args.size() == 3) target = Float3{ std::stof(args[0]), std::stof(args[1]), std::stof(args[2]) };
                else target = Float3{ std::stof(args[0]) };
            }
            catch (...) { throw CliError(std::format("option '{}' contains an invalid number", name)); }
        }, allowInitFrom1 ? 1 : 3, 3);
    }

    void AddOption(const std::vector<std::string>& aliases, bool required, std::vector<int>& target) {
        AddOption(aliases, required, [&target, name = aliases[0]](const auto& args) {
            try { for (const auto& arg : args) target.push_back(std::stoi(arg)); }
            catch (...) { throw CliError(std::format("option '{}' contains an invalid integer", name)); }
        });
    }

    void AddFlag(const std::vector<std::string>& aliases, FlagCallback callback) {
        ValidateAliases(aliases);
        for (const auto& alias : aliases) aliasMap[alias] = aliases[0];
        flags[aliases[0]] = std::move(callback);
    }

    void Parse(const std::vector<std::string>& inputArgs) {
        std::unordered_map<std::string, bool> seen;
        bool optionsEnded = false;
        for (std::size_t i = 2; i < inputArgs.size(); ++i) {
            std::string key = inputArgs[i];
            if (key == "--") { optionsEnded = true; continue; }
            if (optionsEnded) throw CliError(std::format("unexpected positional argument '{}'", key));

            std::optional<std::string> attachedValue;
            if (const auto equals = key.find('='); equals != std::string::npos) {
                attachedValue = key.substr(equals + 1);
                key.resize(equals);
            }

            const auto alias = aliasMap.find(key);
            if (alias == aliasMap.end()) throw CliError(std::format("unrecognized option '{}'", key));
            const std::string& name = alias->second;
            if (name == "--help") {
                if (attachedValue) throw CliError("option '--help' does not take a value");
                throw HelpRequested(helpText);
            }

            if (const auto option = options.find(name); option != options.end()) {
                seen[name] = true;
                std::vector<std::string> args;
                if (attachedValue) args.push_back(std::move(*attachedValue));
                while (args.size() < option->second.maxValues && i + 1 < inputArgs.size()) {
                    const std::string& candidate = inputArgs[i + 1];
                    if (candidate == "--" || IsKnownOption(candidate) || LooksLikeUnknownOption(candidate)) break;
                    args.push_back(candidate);
                    ++i;
                }
                if (args.size() < option->second.minValues)
                    throw CliError(std::format("option '{}' expects a value", name));
                option->second.callback(args);
            }
            else {
                if (attachedValue) throw CliError(std::format("flag '{}' does not take a value", name));
                flags.at(name)();
            }
        }

        for (const auto& [name, option] : options)
            if (option.required && !seen[name]) throw CliError(std::format("missing required option '{}'", name));
    }

    void Parse(int argc, char** argv) { Parse(std::vector<std::string>(argv, argv + argc)); }

private:
    struct Option { bool required; std::size_t minValues; std::size_t maxValues; Callback callback; };

    static void ValidateAliases(const std::vector<std::string>& aliases) {
        if (aliases.empty()) throw std::logic_error("a CLI option must have at least one name");
    }

    bool IsKnownOption(const std::string& candidate) const {
        const auto equals = candidate.find('=');
        return aliasMap.contains(candidate.substr(0, equals));
    }

    static bool LooksLikeUnknownOption(const std::string& candidate) {
        if (candidate.size() < 2 || candidate[0] != '-') return false;
        const char next = candidate[1];
        return next == '-' || ((next < '0' || next > '9') && next != '.');
    }

    std::unordered_map<std::string, std::string> aliasMap;
    std::unordered_map<std::string, Option> options;
    std::unordered_map<std::string, FlagCallback> flags;
    std::string helpText;
};
