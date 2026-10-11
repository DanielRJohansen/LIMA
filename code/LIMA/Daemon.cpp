#include "Daemon.h"

#include "Environment.h"
#include "Format.h"

#include <array>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <streambuf>
#include <thread>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#define WIN32_LEAN_AND_MEAN
#include <winsock2.h>
#include <afunix.h>
#include <windows.h>
#else
#include <fcntl.h>
#include <signal.h>
#include <sys/file.h>
#include <sys/select.h>
#include <sys/socket.h>
#include <sys/stat.h>
#include <sys/un.h>
#include <sys/wait.h>
#include <unistd.h>
#endif

namespace fs = std::filesystem;
using namespace std::chrono_literals;
using Clock = std::chrono::steady_clock;

namespace {

// ------------------------------------------------ Platform ------------------------------------------------ //
// Unix domain sockets exist on both Linux and Windows (10 1803+), so only process spawning, file locking and
// the socket handle type differ between the two.

#ifdef _WIN32
using NativeSocket = SOCKET;
const NativeSocket invalidSocket = INVALID_SOCKET;
void CloseNative(NativeSocket s) { closesocket(s); }
void InitSockets() {
	static const bool initialized = [] { WSADATA data; return WSAStartup(MAKEWORD(2, 2), &data) == 0; }();
	if (!initialized) throw std::runtime_error("Failed to initialize Windows sockets");
}
constexpr int sendFlags = 0;
#else
using NativeSocket = int;
const NativeSocket invalidSocket = -1;
void CloseNative(NativeSocket s) { ::close(s); }
void InitSockets() { signal(SIGPIPE, SIG_IGN); }	// A vanished client must not kill the daemon
constexpr int sendFlags = MSG_NOSIGNAL;
#endif

class Socket {
public:
	Socket() = default;
	explicit Socket(NativeSocket handle) : handle(handle) {}
	Socket(Socket&& other) noexcept : handle(std::exchange(other.handle, invalidSocket)) {}
	Socket& operator=(Socket&& other) noexcept { Close(); handle = std::exchange(other.handle, invalidSocket); return *this; }
	~Socket() { Close(); }

	bool Valid() const { return handle != invalidSocket; }
	NativeSocket Native() const { return handle; }
	void Close() { if (Valid()) CloseNative(std::exchange(handle, invalidSocket)); }

	bool SendAll(const char* data, size_t size) {
		while (size > 0) {
			const int sent = ::send(handle, data, static_cast<int>(std::min<size_t>(size, 1 << 20)), sendFlags);
			if (sent <= 0) return false;
			data += sent;
			size -= sent;
		}
		return true;
	}
	bool ReceiveAll(char* data, size_t size) {
		while (size > 0) {
			const int received = ::recv(handle, data, static_cast<int>(std::min<size_t>(size, 1 << 20)), 0);
			if (received <= 0) return false;
			data += received;
			size -= received;
		}
		return true;
	}
	// Whether data or a connection is waiting, within the timeout
	bool WaitReadable(std::chrono::milliseconds timeout) const {
		fd_set set;
		FD_ZERO(&set);
		FD_SET(handle, &set);
		timeval tv{ static_cast<long>(timeout.count() / 1000), static_cast<long>(timeout.count() % 1000 * 1000) };
		return ::select(static_cast<int>(handle) + 1, &set, nullptr, nullptr, &tv) > 0;
	}

private:
	NativeSocket handle = invalidSocket;
};

sockaddr_un MakeAddress(const fs::path& path) {
	sockaddr_un address{};
	address.sun_family = AF_UNIX;
	const std::string name = path.string();
	if (name.size() >= sizeof(address.sun_path))
		throw std::runtime_error(Lima::Format("The daemon socket path is too long for a Unix socket: {}", name));
	std::memcpy(address.sun_path, name.c_str(), name.size() + 1);
	return address;
}

std::optional<Socket> TryConnect(const fs::path& path) {
	InitSockets();
	Socket socket{ ::socket(AF_UNIX, SOCK_STREAM, 0) };
	if (!socket.Valid()) return std::nullopt;
	const sockaddr_un address = MakeAddress(path);
	if (::connect(socket.Native(), reinterpret_cast<const sockaddr*>(&address), sizeof(address)) != 0)
		return std::nullopt;
	return socket;
}

Socket Listen(const fs::path& path) {
	InitSockets();
	std::error_code ignored;
	fs::remove(path, ignored);	// Left behind by a daemon that did not exit cleanly
	Socket socket{ ::socket(AF_UNIX, SOCK_STREAM, 0) };
	const sockaddr_un address = MakeAddress(path);
	if (!socket.Valid()
		|| ::bind(socket.Native(), reinterpret_cast<const sockaddr*>(&address), sizeof(address)) != 0
		|| ::listen(socket.Native(), SOMAXCONN) != 0)
		throw std::runtime_error(Lima::Format("Failed to listen on {}", path.string()));
	return socket;
}

fs::path ExecutablePath() {
#ifdef _WIN32
	std::wstring buffer(32768, L'\0');
	const DWORD length = GetModuleFileNameW(nullptr, buffer.data(), static_cast<DWORD>(buffer.size()));
	buffer.resize(length);
	return buffer;
#else
	return fs::read_symlink("/proc/self/exe");
#endif
}

int CurrentProcessId() {
#ifdef _WIN32
	return static_cast<int>(GetCurrentProcessId());
#else
	return static_cast<int>(getpid());
#endif
}

// Holds an exclusive lock on a file, so only one client at a time starts a daemon
class FileLock {
public:
	explicit FileLock(const fs::path& path) {
#ifdef _WIN32
		// Share mode 0 makes the open itself the lock
		const auto deadline = Clock::now() + 60s;
		while ((handle = CreateFileW(path.wstring().c_str(), GENERIC_READ | GENERIC_WRITE, 0, nullptr, OPEN_ALWAYS,
			FILE_ATTRIBUTE_NORMAL, nullptr)) == INVALID_HANDLE_VALUE) {
			if (GetLastError() != ERROR_SHARING_VIOLATION || Clock::now() > deadline)
				throw std::runtime_error(Lima::Format("Failed to lock {}", path.string()));
			std::this_thread::sleep_for(20ms);
		}
#else
		fd = ::open(path.c_str(), O_RDWR | O_CREAT | O_CLOEXEC, 0600);
		if (fd < 0 || ::flock(fd, LOCK_EX) != 0)
			throw std::runtime_error(Lima::Format("Failed to lock {}", path.string()));
#endif
	}
	FileLock(const FileLock&) = delete;
	FileLock& operator=(const FileLock&) = delete;
	~FileLock() {
#ifdef _WIN32
		CloseHandle(handle);
#else
		::close(fd);	// Releases the flock
#endif
	}
private:
#ifdef _WIN32
	HANDLE handle = INVALID_HANDLE_VALUE;
#else
	int fd = -1;
#endif
};

// Starts 'lima __daemon SOCKET IDLE' detached from this process and its console, inheriting no handles, so
// it outlives the client and does not hold the client's output pipes open. Returns whether it is still running,
// where that can be known
std::function<bool()> SpawnDetached(const fs::path& executable, const fs::path& socketPath, int idleSeconds, const fs::path& workDir) {
	const std::string idle = std::to_string(idleSeconds);
#ifdef _WIN32
	std::wstring commandline = L"\"" + executable.wstring() + L"\" " + std::wstring(Cli::Daemon::internalServeCommand.begin(),
		Cli::Daemon::internalServeCommand.end()) + L" \"" + socketPath.wstring() + L"\" " + std::wstring(idle.begin(), idle.end());
	STARTUPINFOW startupInfo{ sizeof(startupInfo) };
	PROCESS_INFORMATION processInfo{};
	const DWORD flags = DETACHED_PROCESS | CREATE_NEW_PROCESS_GROUP;
	// Leave the client's job if allowed, so eg. closing an IDE terminal does not kill a running simulation.
	// Jobs that forbid it (like limaclitest's) keep the daemon, which then dies with the job
	BOOL created = CreateProcessW(executable.wstring().c_str(), commandline.data(), nullptr, nullptr, FALSE,
		flags | CREATE_BREAKAWAY_FROM_JOB, nullptr, workDir.wstring().c_str(), &startupInfo, &processInfo);
	if (!created)
		created = CreateProcessW(executable.wstring().c_str(), commandline.data(), nullptr, nullptr, FALSE,
			flags, nullptr, workDir.wstring().c_str(), &startupInfo, &processInfo);
	if (!created)
		throw std::runtime_error(Lima::Format("Failed to start the lima daemon (error {})", GetLastError()));
	CloseHandle(processInfo.hThread);
	const std::shared_ptr<void> process{ processInfo.hProcess, CloseHandle };
	return [process] { return WaitForSingleObject(process.get(), 0) == WAIT_TIMEOUT; };
#else
	const pid_t child = fork();
	if (child < 0) throw std::runtime_error("Failed to start the lima daemon (fork failed)");
	if (child == 0) {
		// Double fork, so the daemon is not our child and is never left as a zombie
		setsid();
		if (fork() != 0) _exit(0);
		const int null = ::open("/dev/null", O_RDWR);
		dup2(null, 0); dup2(null, 1); dup2(null, 2);
		for (int fd = 3; fd < 1024; ++fd) ::close(fd);
		if (chdir(workDir.c_str()) != 0) _exit(127);
		const std::string exe = executable.string(), socket = socketPath.string();
		const std::string command{ Cli::Daemon::internalServeCommand };
		execl(exe.c_str(), exe.c_str(), command.c_str(), socket.c_str(), idle.c_str(), static_cast<char*>(nullptr));
		_exit(127);
	}
	waitpid(child, nullptr, 0);
	return [] { return true; };	// The daemon is our grandchild, so we cannot wait on it
#endif
}

bool WaitForProcessExit(int pid, std::chrono::seconds timeout) {
#ifdef _WIN32
	HANDLE process = OpenProcess(SYNCHRONIZE, FALSE, static_cast<DWORD>(pid));
	if (!process) return true;
	const bool exited = WaitForSingleObject(process, static_cast<DWORD>(timeout.count() * 1000)) == WAIT_OBJECT_0;
	CloseHandle(process);
	return exited;
#else
	const auto deadline = Clock::now() + timeout;
	while (kill(pid, 0) == 0) {
		if (Clock::now() > deadline) return false;
		std::this_thread::sleep_for(50ms);
	}
	return true;
#endif
}

std::string GetEnv(const char* name) {
	const char* value = std::getenv(name);
	return value ? value : "";
}

// ------------------------------------------------ Daemon identity ------------------------------------------------ //

// Each build of lima gets its own daemon, so a rebuilt client never talks to an outdated daemon.
// An old daemon simply idles out
std::string BuildId() {
	const fs::path executable = ExecutablePath();
	const std::string identity = Lima::Format("{}|{}|{}", executable.string(),
		fs::last_write_time(executable).time_since_epoch().count(), fs::file_size(executable));
	uint64_t hash = 14695981039346656037ull;	// FNV-1a
	for (const unsigned char c : identity) { hash ^= c; hash *= 1099511628211ull; }
	return Lima::Format("{:016x}", hash);
}

fs::path RuntimeDir() {
	fs::path dir;
	if (const std::string custom = GetEnv("LIMA_DAEMON_DIR"); !custom.empty())
		dir = custom;
	else {
#ifdef _WIN32
		// Not %LOCALAPPDATA%: packaged (MSIX) apps get it virtualized, which breaks Unix socket files there
		dir = fs::temp_directory_path() / "lima";
#else
		const std::string runtime = GetEnv("XDG_RUNTIME_DIR");
		dir = runtime.empty() ? fs::path(Lima::Format("/tmp/lima-{}", getuid())) : fs::path(runtime) / "lima";
#endif
	}
	fs::create_directories(dir);
#ifndef _WIN32
	fs::permissions(dir, fs::perms::owner_all, fs::perm_options::replace);	// Other users must not submit jobs as us
#endif
	return fs::absolute(dir);
}

struct DaemonPaths {
	fs::path socket, lock, log;
};

DaemonPaths PathsFor(const fs::path& socket) {
	return { socket, fs::path(socket).replace_extension(".lock"), fs::path(socket).replace_extension(".log") };
}

DaemonPaths CurrentPaths() {
	return PathsFor(RuntimeDir() / ("lima-" + BuildId() + ".sock"));
}

int IdleSeconds() {
	const std::string value = GetEnv("LIMA_DAEMON_IDLE");
	try { return value.empty() ? 120 : std::max(1, std::stoi(value)); }
	catch (...) { return 120; }
}

// ------------------------------------------------ Protocol ------------------------------------------------ //
// Request:  a list of strings: protocol version, kind (run/status/stop), the client's working directory, arguments.
// Response: a stream of frames, each a type byte, a 32-bit length and a payload:
//   'A' accepted: from here on the command counts as started, so the client must not retry it
//   'O' / 'E' text for stdout / stderr, 'P' the daemon's process id, 'X' the exit code; always last

constexpr std::string_view protocolVersion = "1";
constexpr uint32_t maxMessageSize = 64u << 20;

void AppendU32(std::string& buffer, uint32_t value) {
	for (int i = 0; i < 4; ++i) buffer.push_back(static_cast<char>((value >> (8 * i)) & 0xFF));
}

bool ReceiveU32(Socket& socket, uint32_t& value) {
	unsigned char bytes[4];
	if (!socket.ReceiveAll(reinterpret_cast<char*>(bytes), 4)) return false;
	value = bytes[0] | bytes[1] << 8 | bytes[2] << 16 | static_cast<uint32_t>(bytes[3]) << 24;
	return true;
}

bool SendStrings(Socket& socket, const std::vector<std::string>& strings) {
	std::string buffer;
	AppendU32(buffer, static_cast<uint32_t>(strings.size()));
	for (const auto& s : strings) {
		AppendU32(buffer, static_cast<uint32_t>(s.size()));
		buffer += s;
	}
	return socket.SendAll(buffer.data(), buffer.size());
}

bool ReceiveStrings(Socket& socket, std::vector<std::string>& strings) {
	uint32_t count = 0;
	if (!ReceiveU32(socket, count) || count > 65536) return false;
	size_t total = 0;
	strings.resize(count);
	for (auto& s : strings) {
		uint32_t size = 0;
		if (!ReceiveU32(socket, size) || (total += size) > maxMessageSize) return false;
		s.resize(size);
		if (!socket.ReceiveAll(s.data(), size)) return false;
	}
	return true;
}

bool SendFrame(Socket& socket, char type, std::string_view payload) {
	std::string buffer(1, type);
	AppendU32(buffer, static_cast<uint32_t>(payload.size()));
	buffer += payload;
	return socket.SendAll(buffer.data(), buffer.size());
}

bool ReceiveFrame(Socket& socket, char& type, std::string& payload) {
	uint32_t size = 0;
	if (!socket.ReceiveAll(&type, 1) || !ReceiveU32(socket, size) || size > maxMessageSize) return false;
	payload.resize(size);
	return socket.ReceiveAll(payload.data(), size);
}

// ------------------------------------------------ Client ------------------------------------------------ //

// Connects to the running daemon, or starts one. nullopt if no daemon could be reached
std::optional<Socket> ConnectOrSpawn(const DaemonPaths& paths, std::string_view fallback) {
	if (auto socket = TryConnect(paths.socket)) return socket;
	FileLock lock{ paths.lock };
	if (auto socket = TryConnect(paths.socket)) return socket;	// Another client started it while we waited

	const auto daemonRunning = SpawnDetached(ExecutablePath(), paths.socket, IdleSeconds(), paths.socket.parent_path());
	// The daemon listens before it initializes anything else, so this is only process startup
	const auto deadline = Clock::now() + 10s;
	while (Clock::now() < deadline && daemonRunning()) {
		if (auto socket = TryConnect(paths.socket)) return socket;
		std::this_thread::sleep_for(25ms);
	}
	std::cerr << "lima: the lima daemon did not start, see " << paths.log.string() << fallback << "\n";
	return std::nullopt;
}

enum class Delivery { NotAccepted, Completed, LostAfterAccepted };

// Sends one request and relays the response. exitCode and pid are set from the response
Delivery Exchange(Socket& socket, const std::vector<std::string>& request, std::ostream* out, std::ostream* err,
	int& exitCode, std::optional<int>* pid = nullptr) {
	if (!SendStrings(socket, request)) return Delivery::NotAccepted;
	bool accepted = false;
	char type = 0;
	std::string payload;
	while (ReceiveFrame(socket, type, payload)) {
		switch (type) {
		case 'A': accepted = true; break;
		case 'O': if (out) { out->write(payload.data(), payload.size()); out->flush(); } break;
		case 'E': if (err) { err->write(payload.data(), payload.size()); err->flush(); } break;
		case 'P': if (pid) *pid = std::stoi(payload); break;
		case 'X': exitCode = std::stoi(payload); return Delivery::Completed;
		default: break;
		}
	}
	return accepted ? Delivery::LostAfterAccepted : Delivery::NotAccepted;
}

std::vector<std::string> MakeRequest(std::string_view kind, const fs::path& cwd, const std::vector<std::string>& args) {
	std::vector<std::string> request{ std::string(protocolVersion), std::string(kind), cwd.string() };
	request.insert(request.end(), args.begin(), args.end());
	return request;
}

// ------------------------------------------------ Server ------------------------------------------------ //

// One client connection. Shared by the stdout and stderr streams of the command it runs
struct Connection {
	Socket socket;
	std::mutex mutex;
	bool alive = true;	// Set when the client has gone; the command still runs to completion

	void Send(char type, std::string_view payload) {
		const std::lock_guard lock(mutex);
		if (alive && !SendFrame(socket, type, payload)) alive = false;
	}
};

// Sends everything written to it to the client, line by line so progress shows up live
class FrameStreambuf : public std::streambuf {
public:
	FrameStreambuf(Connection& connection, char type) : connection(connection), type(type) {}
	~FrameStreambuf() override { sync(); }

protected:
	int_type overflow(int_type c) override {
		if (c != traits_type::eof()) {
			buffer.push_back(traits_type::to_char_type(c));
			if (c == '\n' || buffer.size() >= 4096) sync();
		}
		return traits_type::not_eof(c);
	}
	std::streamsize xsputn(const char* data, std::streamsize size) override {
		buffer.append(data, static_cast<size_t>(size));
		if (std::string_view(data, static_cast<size_t>(size)).find('\n') != std::string_view::npos || buffer.size() >= 4096)
			sync();
		return size;
	}
	int sync() override {
		if (!buffer.empty()) {
			connection.Send(type, buffer);
			buffer.clear();
		}
		return 0;
	}

private:
	Connection& connection;
	char type;
	std::string buffer;
};

std::string Timestamp() {
	const std::time_t now = std::time(nullptr);
	std::tm local{};
#ifdef _WIN32
	localtime_s(&local, &now);
#else
	localtime_r(&now, &local);
#endif
	char text[32];
	std::strftime(text, sizeof(text), "%Y-%m-%d %H:%M:%S", &local);
	return text;
}

// The command line as the user typed it; args[0] is the client's full executable path
std::string JoinArgs(const std::vector<std::string>& args) {
	std::string joined = "lima";
	for (size_t i = 1; i < args.size(); ++i) joined += " " + args[i];
	return joined;
}

class Server {
public:
	Server(fs::path socketPath, std::chrono::seconds idleTimeout, const Cli::Daemon::CommandRunner& runner)
		: socketPath(std::move(socketPath)), idleTimeout(idleTimeout), runner(runner) {}

	int Run() {
		Socket listener = Listen(socketPath);
		Log(Lima::Format("listening on {}, pid {}, idle timeout {} s", socketPath.string(), CurrentProcessId(), idleTimeout.count()));
		Environment::Get();	// Starts creating the CUDA context while the first client connects

		while (true) {
			if (listener.WaitReadable(250ms)) {
				Socket client{ ::accept(listener.Native(), nullptr, nullptr) };
				if (client.Valid()) {
					const std::lock_guard lock(mutex);
					++active;	// Counted before the next idle check, so an accepted client is never dropped
					std::thread([this, client = std::move(client)]() mutable { Handle(std::move(client)); }).detach();
				}
			}
			const std::lock_guard lock(mutex);
			if (stopRequested) { Log("stopping on request"); break; }
			// Results nobody has collected with 'lima wait' keep the daemon alive longer, so a late wait still
			// sees a failure
			const bool uncollected = std::ranges::any_of(jobs, [](const auto& entry) { return entry.second.done; });
			if (active == 0 && Clock::now() - lastActivity >= (uncollected ? std::max(idleTimeout, std::chrono::seconds(3600)) : idleTimeout)) {
				Log("stopping after idle timeout");
				break;
			}
		}

		// Unlink before closing: once closed, a client may start a new daemon that binds the same path, and a
		// later remove would delete its socket. Clients still in the listen backlog get a closed connection
		// before acceptance, and start a new daemon
		std::error_code ignored;
		fs::remove(socketPath, ignored);
		listener.Close();

		std::unique_lock lock(mutex);
		if (active > 0) Log(Lima::Format("waiting for {} running commands", active));
		changed.wait(lock, [this] { return active == 0; });
		Log("exited");
		return 0;
	}

private:
	struct Job {
		std::string command;
		bool async = false;
		bool done = false;
		int exitCode = 0;
		fs::path log;	// Output of an async job
	};

	void Handle(Socket socket) {
		uint64_t id = 0;
		try {
			Connection connection{ std::move(socket) };
			std::vector<std::string> request;
			if (!ReceiveStrings(connection.socket, request) || request.size() < 3) {
				Finish();
				return;
			}
			if (request[0] != protocolVersion) {
				connection.Send('E', "lima: the client and the daemon are different versions\n");
				connection.Send('X', "1");
				Finish();
				return;
			}
			const std::string kind = request[1];
			const fs::path cwd = request[2];
			const std::vector<std::string> args(request.begin() + 3, request.end());
			connection.Send('A', "");

			if (kind == "run") {
				id = AddJob(args, false, {});
				Log(Lima::Format("[{}] start in {}: {}", id, cwd.string(), JoinArgs(args)));
				// A client that goes away (eg. Ctrl+C) does not stop its command; it still completes and writes its output
				int exitCode = 1;
				{
					FrameStreambuf outBuffer{ connection, 'O' }, errBuffer{ connection, 'E' };
					std::ostream out{ &outBuffer }, err{ &errBuffer };
					exitCode = runner(Cli::Context{ cwd, out, err }, args);
				}
				Log(Lima::Format("[{}] exit {}{}", id, exitCode, connection.alive ? "" : " (client disconnected)"));
				RemoveJob(id);
				connection.Send('X', std::to_string(exitCode));
			}
			else if (kind == "submit") {
				const uint64_t jobId = NextId();
				const fs::path log = cwd / Lima::Format("lima-{}.log", jobId);
				AddJob(args, true, log, jobId);
				Log(Lima::Format("[{}] submitted in {}: {}", jobId, cwd.string(), JoinArgs(args)));
				connection.Send('O', Lima::Format("Submitted job {}, output in {}\n", jobId, log.string()));
				connection.Send('X', "0");
				connection.socket.Close();	// The client returns now

				int exitCode = 1;
				std::ofstream file(log);
				if (file) {
					exitCode = runner(Cli::Context{ cwd, file, file }, args);
					file << Lima::Format("\nlima: exit code {}\n", exitCode);
				}
				else
					Log(Lima::Format("[{}] cannot write {}", jobId, log.string()));
				Log(Lima::Format("[{}] exit {}", jobId, exitCode));
				CompleteJob(jobId, exitCode);
			}
			else if (kind == "wait") {
				const auto [text, exitCode] = WaitForAsyncJobs();
				connection.Send('O', text);
				connection.Send('X', std::to_string(exitCode));
			}
			else if (kind == "status") {
				connection.Send('O', StatusText());
				connection.Send('X', "0");
			}
			else if (kind == "stop") {
				connection.Send('P', std::to_string(CurrentProcessId()));
				connection.Send('O', StatusText());
				{
					const std::lock_guard lock(mutex);
					stopRequested = true;
				}
				connection.Send('X', "0");
			}
			else {
				connection.Send('E', Lima::Format("lima: unknown daemon request '{}'\n", kind));
				connection.Send('X', "1");
			}
		}
		catch (const std::exception& error) {
			Log(Lima::Format("[{}] error: {}", id, error.what()));
		}
		Finish();
	}

	uint64_t NextId() {
		const std::lock_guard lock(mutex);
		return ++nextId;
	}

	uint64_t AddJob(const std::vector<std::string>& args, bool async, fs::path log, uint64_t id = 0) {
		const std::lock_guard lock(mutex);
		if (id == 0) id = ++nextId;
		jobs[id] = Job{ JoinArgs(args), async, false, 0, std::move(log) };
		return id;
	}

	void RemoveJob(uint64_t id) {
		const std::lock_guard lock(mutex);
		jobs.erase(id);
	}

	void CompleteJob(uint64_t id, int exitCode) {
		{
			const std::lock_guard lock(mutex);
			jobs[id].done = true;
			jobs[id].exitCode = exitCode;
		}
		changed.notify_all();
	}

	// Waits for every async job, then reports and forgets the finished ones. Exit code 1 if any failed
	std::pair<std::string, int> WaitForAsyncJobs() {
		std::unique_lock lock(mutex);
		changed.wait(lock, [this] {
			return std::ranges::all_of(jobs, [](const auto& entry) { return !entry.second.async || entry.second.done; });
		});
		std::string text;
		int failed = 0;
		for (auto it = jobs.begin(); it != jobs.end();) {
			const Job& job = it->second;
			if (!job.async) { ++it; continue; }
			text += Lima::Format("[{}] {}  {}\n", it->first, job.exitCode == 0 ? "done  " : Lima::Format("FAILED (exit {})", job.exitCode), job.command);
			if (job.exitCode != 0) {
				++failed;
				text += Lima::Format("      see {}\n", job.log.string());
			}
			it = jobs.erase(it);
		}
		if (text.empty()) text = "No jobs to wait for\n";
		else if (failed > 0) text += Lima::Format("{} job{} failed\n", failed, failed == 1 ? "" : "s");
		return { text, failed > 0 ? 1 : 0 };
	}

	void Finish() {
		{
			const std::lock_guard lock(mutex);
			--active;
			lastActivity = Clock::now();
		}
		changed.notify_all();
	}

	std::string StatusText() {
		const std::lock_guard lock(mutex);
		const auto running = std::ranges::count_if(jobs, [](const auto& entry) { return !entry.second.done; });
		std::string text = Lima::Format("lima daemon, pid {}, {} running command{}\n", CurrentProcessId(), running,
			running == 1 ? "" : "s");
		for (const auto& [id, job] : jobs)
			text += Lima::Format("  [{}] {}{}\n", id, job.done ? Lima::Format("finished (exit {}), not yet collected by 'lima wait': ", job.exitCode) : "", job.command);
		return text;
	}

	void Log(const std::string& message) {
		const std::lock_guard lock(logMutex);
		std::cout << Timestamp() << "  " << message << std::endl;
	}

	fs::path socketPath;
	std::chrono::seconds idleTimeout;
	const Cli::Daemon::CommandRunner& runner;

	std::mutex mutex;
	std::condition_variable changed;	// A connection or a job finished
	size_t active = 0;
	uint64_t nextId = 0;
	std::map<uint64_t, Job> jobs;	// Running commands, and async results not yet collected by 'lima wait'
	Clock::time_point lastActivity = Clock::now();
	bool stopRequested = false;
	std::mutex logMutex;
};

bool IsAnyOf(std::string_view arg, std::initializer_list<std::string_view> options) {
	return std::ranges::find(options, arg) != options.end();
}

} // namespace

// ------------------------------------------------ Public ------------------------------------------------ //

bool Cli::Daemon::ShouldUse(const std::vector<std::string>& args) {
	const std::string noDaemon = GetEnv("LIMA_NO_DAEMON");
	if (!noDaemon.empty() && noDaemon != "0") return false;
	// Only GPU commands gain from batching. Displays stay in the user's own process
	if (args.size() < 2 || !IsAnyOf(args[1], { "mdrun", "em", "buildmembrane" })) return false;
	return std::ranges::none_of(args.begin() + 2, args.end(), [](const std::string& arg) {
		return IsAnyOf(arg, { "-d", "--display", "-display", "-h", "--help", "-help" });
	});
}

std::optional<int> Cli::Daemon::RunRemote(const std::vector<std::string>& args, bool async) {
	std::vector<std::string> request;
	try {
		request = MakeRequest(async ? "submit" : "run", fs::current_path(), args);
	}
	catch (const std::exception&) {
		// Eg. a working directory outside the system code page. Running in-process reports the real error
		return std::nullopt;
	}
	const DaemonPaths paths = CurrentPaths();
	const std::string_view fallback = async ? "" : ", running in this process instead";
	// A second attempt covers a daemon that was shutting down after its idle timeout as we connected
	for (int attempt = 0; attempt < 2; ++attempt) {
		auto socket = ConnectOrSpawn(paths, fallback);
		if (!socket) return std::nullopt;
		int exitCode = 1;
		switch (Exchange(*socket, request, &std::cout, &std::cerr, exitCode)) {
		case Delivery::Completed: return exitCode;
		case Delivery::LostAfterAccepted:
			std::cerr << "lima: lost the connection to the lima daemon, see " << paths.log.string() << "\n";
			return 1;
		case Delivery::NotAccepted: break;
		}
	}
	std::cerr << "lima: the lima daemon did not accept the command" << fallback << "\n";
	return std::nullopt;
}

int Cli::Daemon::RunWait(const Context& ctx, const std::vector<std::string>&) {
	auto socket = TryConnect(CurrentPaths().socket);
	if (!socket) {
		ctx.out << "No lima daemon is running, so there are no jobs to wait for\n";
		return 0;
	}
	int exitCode = 1;
	if (Exchange(*socket, MakeRequest("wait", ctx.cwd, {}), &ctx.out, &ctx.err, exitCode) != Delivery::Completed) {
		ctx.err << "lima: lost the connection to the lima daemon\n";
		return 1;
	}
	return exitCode;
}

int Cli::Daemon::RunDaemonCommand(const Context& ctx, const std::vector<std::string>& args) {
	if (args.size() != 3 || !IsAnyOf(args[2], { "status", "stop" }))
		throw std::invalid_argument("usage: lima daemon status|stop");
	const DaemonPaths paths = CurrentPaths();
	auto socket = TryConnect(paths.socket);
	if (!socket) {
		ctx.out << "No lima daemon is running\n";
		return 0;
	}
	int exitCode = 1;
	std::optional<int> pid;
	if (Exchange(*socket, MakeRequest(args[2], ctx.cwd, {}), &ctx.out, &ctx.err, exitCode, &pid) != Delivery::Completed) {
		ctx.err << "lima: lost the connection to the lima daemon\n";
		return 1;
	}
	if (args[2] == "stop" && pid) {
		ctx.out << "Waiting for the daemon to finish its running commands and exit\n" << std::flush;
		if (!WaitForProcessExit(*pid, 3600s)) {
			ctx.err << "lima: the daemon did not exit\n";
			return 1;
		}
		ctx.out << "Stopped\n";
	}
	return exitCode;
}

int Cli::Daemon::Serve(const std::vector<std::string>& args, const CommandRunner& runner) {
	if (args.size() != 4) throw std::invalid_argument("usage: lima __daemon SOCKET IDLE_SECONDS");
	const fs::path socketPath = args[2];
	const DaemonPaths paths = PathsFor(socketPath);
	// A detached process has no console. Everything the daemon itself prints goes to its log. Appended, so a
	// replacement daemon keeps the log of one that crashed
	if (!std::freopen(paths.log.string().c_str(), "a", stdout) || !std::freopen(paths.log.string().c_str(), "a", stderr))
		throw std::runtime_error("Failed to open the daemon log");
	return Server{ socketPath, std::chrono::seconds(std::stoi(args[3])), runner }.Run();
}
