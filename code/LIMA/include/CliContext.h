#pragma once

#include <filesystem>
#include <ostream>

namespace Cli {

// Everything a command may take from its process. Commands run either in their own process or on a thread of the
// lima daemon, where several run at once, so they must use this instead of fs::current_path(), std::cout or printf.
struct Context {
    std::filesystem::path cwd;
    std::ostream& out;
    std::ostream& err;

    // Relative paths are relative to the user's working directory, not the daemon's
    std::filesystem::path Resolve(const std::filesystem::path& path) const {
        return path.empty() || path.is_absolute() ? path : (cwd / path).lexically_normal();
    }
};

}
