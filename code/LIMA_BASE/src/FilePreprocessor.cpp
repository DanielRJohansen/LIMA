
#include <fstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <system_error>
#include <cctype>

#include "Filehandling.h"

namespace fs = std::filesystem;

// ========================= Utilities ============================================================

static inline std::string ReadFile(const fs::path& file) {
    std::error_code errorCode;
    const auto fileSize = fs::file_size(file, errorCode);
    std::ifstream inputStream(file, std::ios::binary);
    if (!inputStream) return {};
    std::string data;
    if (!errorCode && fileSize > 0) data.reserve(static_cast<size_t>(fileSize));
    data.assign(std::istreambuf_iterator<char>(inputStream), std::istreambuf_iterator<char>());
    return data;
}

static inline fs::path MakeWeaklyCanonical(const fs::path& p) {
    std::error_code ec;
    auto canonicalPath = fs::weakly_canonical(p, ec);
    return ec ? fs::absolute(p) : canonicalPath;
}

static inline void TrimInPlace(std::string_view& s) {
    auto isSpace = [](unsigned char c) { return c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == '\f'; };
    while (!s.empty() && isSpace(static_cast<unsigned char>(s.front()))) s.remove_prefix(1);
    while (!s.empty() && isSpace(static_cast<unsigned char>(s.back())))  s.remove_suffix(1);
}

static inline std::string_view TakeUntilEndOfLine(std::string_view& s) {
    size_t pos = s.find('\n');
    std::string_view line = pos == std::string_view::npos ? s : s.substr(0, pos);
    s = pos == std::string_view::npos ? std::string_view{} : s.substr(pos + 1);
    return line;
}

static inline bool LineStartsWithHash(std::string_view s) {
    auto i = s.find_first_not_of(" \t");
    return i != std::string_view::npos && s[i] == '#';
}

static inline std::string JoinLineContinuations(std::string_view text) {
    std::string out;
    out.reserve(text.size());
    for (size_t i = 0; i < text.size(); ++i) {
        char c = text[i];
        if (c == '\\') {
            if (i + 1 < text.size() && text[i + 1] == '\n') { ++i; continue; }
            if (i + 2 < text.size() && text[i + 1] == '\r' && text[i + 2] == '\n') { i += 2; continue; }
        }
        out.push_back(c);
    }
    return out;
}

static inline std::string_view NextIdentifier(std::string_view& s) {
    TrimInPlace(s);
    size_t i = 0;
    while (i < s.size() && (std::isalnum(static_cast<unsigned char>(s[i])) || s[i] == '_')) ++i;
    std::string_view id = s.substr(0, i);
    s.remove_prefix(i);
    return id;
}

static inline fs::path ResolveIncludeTarget(std::string_view spec,
    const fs::path& currentDirectory,
    const std::vector<fs::path>& includeDirectories)
{
    TrimInPlace(spec);
    if (spec.empty()) return {};
    char endMarker = (spec.front() == '"') ? '"' : ((spec.front() == '<') ? '>' : '\0');
    if (!endMarker) return {};
    size_t start = 1;
    size_t end = spec.find(endMarker, start);
    if (end == std::string_view::npos) return {};
    std::string relative(spec.substr(start, end - start));

    // Try local directory first
    fs::path candidate = MakeWeaklyCanonical(currentDirectory / relative);
    if (fs::exists(candidate)) return candidate;

    // Then search include directories
    for (const auto& dir : includeDirectories) {
        fs::path p = MakeWeaklyCanonical(dir / relative);
        if (fs::exists(p)) return p;
    }
    return {};
}

// ========================= Conditional State ====================================================

struct ConditionalFrame {
    bool parentActive;
    bool conditionTrue;
    bool inElse;
};

static inline bool ComputeActive(const std::vector<ConditionalFrame>& stack) {
    if (stack.empty()) return true;
    const auto& f = stack.back();
    const bool thisActive = f.inElse ? (!f.conditionTrue) : f.conditionTrue;
    return f.parentActive && thisActive;
}

// ========================= Core Implementation ==================================================

static std::string PreprocessFileImpl(const fs::path& file,
    const std::vector<fs::path>& includeDirectories,
    std::unordered_set<fs::path>& seen,
    std::unordered_set<std::string>& definedSymbols)
{
    std::string content = ReadFile(file);
    if (content.empty()) return {};

    content = JoinLineContinuations(content);

    std::string_view sv{ content };
    std::string output;
    output.reserve(content.size() + content.size() / 8);

    const fs::path currentDirectory = file.parent_path();
    std::vector<ConditionalFrame> conditionalStack;

    auto lineIsActive = [&]() { return ComputeActive(conditionalStack); };

    while (!sv.empty()) {
        std::string_view rawLine = TakeUntilEndOfLine(sv);
        std::string_view line = rawLine;

        if (!LineStartsWithHash(line)) {
            if (lineIsActive()) {
                output.append(line);
                output.push_back('\n');
            }
            continue;
        }

        // Parse directive
        size_t hashPos = line.find('#');
        std::string_view rest = hashPos == std::string_view::npos ? std::string_view{} : line.substr(hashPos + 1);
        TrimInPlace(rest);
        std::string_view keyword = NextIdentifier(rest);
        TrimInPlace(rest);

        // ----- include -----
        if (keyword == "include") {
            if (lineIsActive()) {
                fs::path includePath = ResolveIncludeTarget(rest, currentDirectory, includeDirectories);
                if (!includePath.empty()) {
                    auto canonicalPath = MakeWeaklyCanonical(includePath);
                    if (!seen.contains(canonicalPath)) {
                        seen.insert(canonicalPath);
                        output += PreprocessFileImpl(canonicalPath, includeDirectories, seen, definedSymbols);
                    }
                }
            }
            continue;
        }

        // ----- define / undefine (boolean symbols only) -----
        if (keyword == "define") {
            if (lineIsActive()) {
                std::string symbol(NextIdentifier(rest));
                if (!symbol.empty()) definedSymbols.insert(std::move(symbol));
            }
            continue;
        }
        if (keyword == "undef" || keyword == "undefine") {
            if (lineIsActive()) {
                std::string symbol(NextIdentifier(rest));
                if (!symbol.empty()) definedSymbols.erase(symbol);
            }
            continue;
        }

        // ----- conditionals -----
        if (keyword == "ifdef") {
            std::string symbol(NextIdentifier(rest));
            bool parentActive = lineIsActive();
            bool conditionTrue = definedSymbols.find(symbol) != definedSymbols.end();
            conditionalStack.push_back({ parentActive, conditionTrue, false });
            continue;
        }
        if (keyword == "ifndef") {
            std::string symbol(NextIdentifier(rest));
            bool parentActive = lineIsActive();
            bool conditionTrue = !(definedSymbols.find(symbol) != definedSymbols.end());
            conditionalStack.push_back({ parentActive, conditionTrue, false });
            continue;
        }
        if (keyword == "else") {
            if (!conditionalStack.empty()) {
                conditionalStack.back().inElse = true;
            }
            continue;
        }
        if (keyword == "endif") {
            if (!conditionalStack.empty()) conditionalStack.pop_back();
            continue;
        }

        // Unknown directive: pass through if active
        if (lineIsActive()) {
            output.push_back('#');
            output.append(keyword);
            if (!rest.empty()) { output.push_back(' '); output.append(rest); }
            output.push_back('\n');
        }
    }

    return output;
}

// ========================= Public API ===========================================================

std::string FileUtils::PreprocessFile(const fs::path& file, const std::vector<fs::path>& includeDirs, std::unordered_set<std::string>& defines) {
    auto canonical = MakeWeaklyCanonical(file);

    std::unordered_set<fs::path> seen{};
    seen.insert(canonical);

    return PreprocessFileImpl(canonical, includeDirs, seen, defines);
}
