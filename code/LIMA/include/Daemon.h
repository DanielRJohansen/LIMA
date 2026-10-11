#pragma once

#include "CliContext.h"

#include <functional>
#include <optional>
#include <string>
#include <vector>

// The lima daemon is a background lima process that owns one Environment. GPU commands are sent to it, so
// commands started at the same time, eg. 'lima mdrun ... &' in a loop, are batched on the GPU together.
// It is started by the first command that needs it and exits after a period without work.
// Only the lima CLI uses it; everything else uses Environment directly.
namespace Cli::Daemon {

// Runs a command line in the calling thread and returns its exit code. Must not throw.
using CommandRunner = std::function<int(const Context&, const std::vector<std::string>& args)>;

// Whether this command line should be sent to the daemon instead of running in this process
bool ShouldUse(const std::vector<std::string>& args);

// Runs the command on the daemon, starting one if none is running, and relays its output. Returns the exit code,
// or nullopt if no daemon could take the command, in which case it has not run and may run in-process instead.
std::optional<int> RunRemote(const std::vector<std::string>& args);

// 'lima daemon status|stop'
int RunDaemonCommand(const Context& ctx, const std::vector<std::string>& args);

// Entry point of the daemon process itself ('lima __daemon SOCKET IDLE_SECONDS')
int Serve(const std::vector<std::string>& args, const CommandRunner& runner);

inline constexpr std::string_view internalServeCommand = "__daemon";

}
