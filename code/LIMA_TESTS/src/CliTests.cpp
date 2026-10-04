// limaclitest: release tests for the lima command line programs.
//
// Every test runs the real lima executable as a child process in a fresh directory, and checks the contract a
// user relies on: exit codes, which files appear, and properties of the output that follow directly from the
// inputs and arguments, with generous tolerances. Exact coordinates, energies and console text are never compared,
// so changes to the simulation or to output formatting should not break these tests.
// Everything that can display does, so the run can be followed visually.
//
// Usage: limaclitest [--lima PATH] [--no-preview] [--dwell SECONDS] [COMMAND...]
//   --lima PATH       lima executable to test (default: the one in this build tree)
//   --no-preview      Dont open 'lima render' on the output of commands that have no display of their own
//   --dwell SECONDS   How long the render test keeps the window open (default: 5)
//   COMMAND...        Only run the tests for these commands, eg. 'limaclitest render mdrun'

#include "TestUtils.h"
#include "Format.h"
#include "CliDefinitions.h"
#include "MDFiles.h"
#include "SimulationBuilder.h"
#include "xdrfile_trr.h"

#include <cctype>
#include <cfloat>
#include <cmath>
#include <set>
#include <unordered_set>

#ifndef NOMINMAX
#define NOMINMAX
#endif
#define WIN32_LEAN_AND_MEAN
#include <windows.h>

#ifndef PW_RENDERFULLCONTENT
#define PW_RENDERFULLCONTENT 0x00000002
#endif

using namespace TestUtils;
using namespace std::chrono_literals;

namespace {

struct Config {
	fs::path lima;
	fs::path fixtures = FileUtils::GetLimaDir() / "tests" / "clitests" / "fixtures";
	fs::path runDir;
	bool preview = true;
	std::chrono::seconds renderDwell{ 5 };
	std::chrono::seconds previewDwell{ 4 };
	std::set<std::string> only;
};
Config config;

// ------------------------------------------------ Child processes ------------------------------------------------ //

// Children are assigned to this job, so they die with the test runner instead of lingering with open windows
HANDLE jobObject = nullptr;

std::wstring ToWide(const std::string& s) {
	if (s.empty()) return {};
	const int n = MultiByteToWideChar(CP_UTF8, 0, s.data(), static_cast<int>(s.size()), nullptr, 0);
	std::wstring out(n, L'\0');
	MultiByteToWideChar(CP_UTF8, 0, s.data(), static_cast<int>(s.size()), out.data(), n);
	return out;
}

std::string Quote(const std::string& arg) {
	return arg.empty() || arg.find_first_of(" \t") != std::string::npos ? "\"" + arg + "\"" : arg;
}

struct RunOptions {
	std::chrono::seconds timeout = 120s;
	// stderr is always captured. stdout is only captured when the test needs it; otherwise it goes straight
	// to the console so progress is shown live
	bool captureStdout = false;
	bool echo = true;	// Show the captured output on the console as it arrives
};

class ChildProcess;
// The child whose windows are being watched. Only one child runs at a time
ChildProcess* windowHookTarget = nullptr;

class ChildProcess {
public:
	ChildProcess(const std::vector<std::string>& args, const fs::path& workDir, bool captureStdout, bool echo = true) {
		std::string commandline = Quote(config.lima.string());
		std::string display = "lima";
		for (const auto& arg : args) {
			commandline += " " + Quote(arg);
			display += " " + Quote(arg);
		}
		if (echo) std::cout << "\033[90m> " << display << "\033[0m\n" << std::flush;

		SECURITY_ATTRIBUTES sa{ sizeof(sa), nullptr, TRUE };
		HANDLE writePipe = nullptr;
		if (!CreatePipe(&readPipe, &writePipe, &sa, 0))
			throw std::runtime_error("CreatePipe failed");
		SetHandleInformation(readPipe, HANDLE_FLAG_INHERIT, 0);

		HANDLE nul = CreateFileW(L"NUL", GENERIC_READ, FILE_SHARE_READ | FILE_SHARE_WRITE, &sa, OPEN_EXISTING, 0, nullptr);
		HANDLE stdoutHandle = nullptr;
		if (!captureStdout)
			DuplicateHandle(GetCurrentProcess(), GetStdHandle(STD_OUTPUT_HANDLE), GetCurrentProcess(), &stdoutHandle, 0, TRUE, DUPLICATE_SAME_ACCESS);

		STARTUPINFOW si{ sizeof(si) };
		si.dwFlags = STARTF_USESTDHANDLES;
		si.hStdInput = nul;
		si.hStdOutput = captureStdout ? writePipe : stdoutHandle;
		si.hStdError = writePipe;

		std::wstring wideCommandline = ToWide(commandline);
		PROCESS_INFORMATION pi{};
		// Suspended, so the window hook is in place before the child can show anything
		const BOOL created = CreateProcessW(nullptr, wideCommandline.data(), nullptr, nullptr, TRUE, CREATE_SUSPENDED, nullptr,
			workDir.wstring().c_str(), &si, &pi);
		const DWORD error = GetLastError();
		CloseHandle(writePipe);	// Only the child may hold the write end, so we get EOF when it exits
		if (nul != INVALID_HANDLE_VALUE) CloseHandle(nul);
		if (stdoutHandle) CloseHandle(stdoutHandle);
		if (!created) {
			CloseHandle(readPipe);
			throw TestFailure{ Lima::Format("Failed to start {} (error {})", config.lima.string(), error) };
		}
		process = pi.hProcess;
		pid = pi.dwProcessId;
		if (jobObject) AssignProcessToJobObject(jobObject, process);
		// Notified of every window the child shows, however briefly. Polling misses windows that live for
		// only a few milliseconds, eg. when a minimization with -d converges immediately
		windowHookTarget = this;
		showHook = SetWinEventHook(EVENT_OBJECT_SHOW, EVENT_OBJECT_SHOW, nullptr, OnWindowShown, pid, 0, WINEVENT_OUTOFCONTEXT);
		ResumeThread(pi.hThread);
		CloseHandle(pi.hThread);
		started = std::chrono::steady_clock::now();

		this->echo = echo;
		reader = std::thread([this, echo] {
			char buffer[4096];
			DWORD n = 0;
			while (ReadFile(readPipe, buffer, sizeof(buffer), &n, nullptr) && n > 0) {
				if (echo) {
					std::fwrite(buffer, 1, n, stdout);
					std::fflush(stdout);
				}
				std::lock_guard lock(outputMutex);
				output.append(buffer, n);
			}
		});
	}

	ChildProcess(const ChildProcess&) = delete;
	ChildProcess& operator=(const ChildProcess&) = delete;

	~ChildProcess() {
		if (showHook) UnhookWinEvent(showHook);
		windowHookTarget = nullptr;
		if (Running()) Kill();
		if (reader.joinable()) reader.join();
		if (echo) std::cout << '\n' << std::flush;	// lima does not always end its output with a newline
		CloseHandle(readPipe);
		CloseHandle(process);
	}

	bool Running() const { return WaitForSingleObject(process, 0) == WAIT_TIMEOUT; }

	// Returns the exit code, or nullopt if the process was still running after the timeout.
	// Pumps this thread's messages while waiting, which is how the window hook is delivered
	std::optional<int> WaitFor(std::chrono::milliseconds timeout) {
		const auto deadline = std::chrono::steady_clock::now() + timeout;
		while (true) {
			PumpMessages();
			const auto remaining = std::max(std::chrono::milliseconds{ 0 },
				std::chrono::duration_cast<std::chrono::milliseconds>(deadline - std::chrono::steady_clock::now()));
			const DWORD wait = MsgWaitForMultipleObjects(1, &process, FALSE, static_cast<DWORD>(remaining.count()), QS_ALLINPUT);
			if (wait == WAIT_OBJECT_0) {
				PumpMessages();
				DWORD code = 0;
				GetExitCodeProcess(process, &code);
				return static_cast<int>(code);
			}
			if (wait != WAIT_OBJECT_0 + 1)	// Timeout or failure. +1 means a message arrived
				return std::nullopt;
		}
	}

	// Whether the child has shown any window so far
	bool WindowShown() const { return windowShown; }

	void Kill() {
		TerminateProcess(process, 0xDEAD);
		WaitForSingleObject(process, 5000);
	}

	// The first visible top-level window owned by the process, if any
	HWND FindOwnWindow() const {
		struct Search { DWORD pid; HWND found; } search{ pid, nullptr };
		EnumWindows([](HWND hwnd, LPARAM param) -> BOOL {
			auto& s = *reinterpret_cast<Search*>(param);
			DWORD owner = 0;
			GetWindowThreadProcessId(hwnd, &owner);
			if (owner == s.pid && IsWindowVisible(hwnd) && GetWindow(hwnd, GW_OWNER) == nullptr) {
				s.found = hwnd;
				return FALSE;
			}
			return TRUE;
		}, reinterpret_cast<LPARAM>(&search));
		return search.found;
	}

	// Captured output; only complete once the process has exited
	std::string Output() {
		if (!Running() && reader.joinable()) reader.join();
		std::lock_guard lock(outputMutex);
		return output;
	}

	std::chrono::duration<double> Elapsed() const { return std::chrono::steady_clock::now() - started; }

private:
	static void CALLBACK OnWindowShown(HWINEVENTHOOK, DWORD, HWND, LONG idObject, LONG idChild, DWORD, DWORD) {
		if (idObject == OBJID_WINDOW && idChild == CHILDID_SELF && windowHookTarget)
			windowHookTarget->windowShown = true;
	}

	static void PumpMessages() {
		MSG message;
		while (PeekMessageW(&message, nullptr, 0, 0, PM_REMOVE)) {
			TranslateMessage(&message);
			DispatchMessageW(&message);
		}
	}

	HWINEVENTHOOK showHook = nullptr;
	bool windowShown = false;
	HANDLE process = nullptr;
	HANDLE readPipe = nullptr;
	DWORD pid = 0;
	std::thread reader;
	bool echo = true;
	std::mutex outputMutex;
	std::string output;
	std::chrono::steady_clock::time_point started;
};

struct ProcessResult {
	std::optional<int> exitCode;	// nullopt if killed after timing out
	bool windowSeen = false;
	std::string output;
	std::string commandline;
};

std::string Tail(const std::string& text, size_t maxChars = 600) {
	std::string tail = text.size() > maxChars ? "..." + text.substr(text.size() - maxChars) : text;
	while (!tail.empty() && std::isspace(static_cast<unsigned char>(tail.back()))) tail.pop_back();
	return tail;
}

// Crashes show up as NTSTATUS codes, which are only recognizable in hex
std::string DescribeExitCode(int code) {
	return code < 0 ? Lima::Format("{} (0x{:08X})", code, static_cast<uint32_t>(code)) : std::to_string(code);
}

ProcessResult Run(const fs::path& workDir, const std::vector<std::string>& args, RunOptions options = {}) {
	ProcessResult result;
	result.commandline = "lima";
	for (const auto& arg : args) result.commandline += " " + Quote(arg);

	ChildProcess child(args, workDir, options.captureStdout, options.echo);
	while (true) {
		if (auto code = child.WaitFor(100ms)) {
			result.exitCode = code;
			break;
		}
		if (child.Elapsed() > options.timeout) {
			child.Kill();
			break;
		}
	}
	result.windowSeen = child.WindowShown();
	result.output = child.Output();
	std::cout << std::flush;
	return result;
}

void RequireExitCode(const ProcessResult& result, int expected) {
	ASSERT(result.exitCode.has_value(), Lima::Format("'{}' timed out and was killed", result.commandline));
	ASSERT(*result.exitCode == expected, Lima::Format("'{}' exited with {}, expected {}. Output:\n{}",
		result.commandline, DescribeExitCode(*result.exitCode), expected, Tail(result.output)));
}

void RequireSuccess(const ProcessResult& result) { RequireExitCode(result, 0); }

void RequireWindowSeen(const ProcessResult& result) {
	ASSERT(result.windowSeen, Lima::Format("'{}' was asked to display, but no window appeared", result.commandline));
}

// ------------------------------------------------ Rendering ------------------------------------------------ //

// Captures the window through DWM (works for occluded OpenGL windows), saves it as a .bmp, and returns the
// number of distinct colors, capped at 256
int CaptureWindow(HWND hwnd, const fs::path& bmpPath) {
	RECT rect{};
	GetWindowRect(hwnd, &rect);
	const int w = rect.right - rect.left;
	const int h = rect.bottom - rect.top;
	if (w <= 0 || h <= 0) return 0;

	HDC screen = GetDC(nullptr);
	HDC memory = CreateCompatibleDC(screen);
	HBITMAP bitmap = CreateCompatibleBitmap(screen, w, h);
	HGDIOBJ previous = SelectObject(memory, bitmap);
	const BOOL printed = PrintWindow(hwnd, memory, PW_RENDERFULLCONTENT);

	BITMAPINFO info{};
	info.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
	info.bmiHeader.biWidth = w;
	info.bmiHeader.biHeight = -h;	// Top-down
	info.bmiHeader.biPlanes = 1;
	info.bmiHeader.biBitCount = 32;
	info.bmiHeader.biCompression = BI_RGB;
	std::vector<uint32_t> pixels(static_cast<size_t>(w) * h);
	GetDIBits(memory, bitmap, 0, h, pixels.data(), &info, DIB_RGB_COLORS);

	SelectObject(memory, previous);
	DeleteObject(bitmap);
	DeleteDC(memory);
	ReleaseDC(nullptr, screen);
	if (!printed) return 0;

	std::ofstream file(bmpPath, std::ios::binary);
	BITMAPFILEHEADER fileHeader{};
	fileHeader.bfType = 0x4D42;	// "BM"
	fileHeader.bfOffBits = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER);
	fileHeader.bfSize = fileHeader.bfOffBits + static_cast<DWORD>(pixels.size() * sizeof(uint32_t));
	file.write(reinterpret_cast<const char*>(&fileHeader), sizeof(fileHeader));
	file.write(reinterpret_cast<const char*>(&info.bmiHeader), sizeof(BITMAPINFOHEADER));
	file.write(reinterpret_cast<const char*>(pixels.data()), pixels.size() * sizeof(uint32_t));

	std::unordered_set<uint32_t> colors;
	for (const uint32_t pixel : pixels) {
		colors.insert(pixel & 0x00FFFFFF);
		if (colors.size() >= 256) break;
	}
	return static_cast<int>(colors.size());
}

struct RenderSession {
	std::string error;	// Empty when the window appeared, stayed alive and responsive, and closed cleanly
	int distinctColors = -1;
	std::chrono::duration<double> timeToWindow{};
};

// lima render runs until its window is closed. Wait for the window, keep it open for a while, then close it
// the way a user would (WM_CLOSE = clicking X) and require a clean exit
RenderSession RunRenderAndClose(const fs::path& workDir, const std::vector<std::string>& args,
	std::chrono::seconds dwell, std::optional<fs::path> screenshot)
{
	RenderSession session;
	ChildProcess child(args, workDir, false);

	HWND window = nullptr;
	while (!(window = child.FindOwnWindow())) {
		if (auto code = child.WaitFor(100ms)) {
			session.error = Lima::Format("render exited with {} before a window appeared. Output:\n{}", DescribeExitCode(*code), Tail(child.Output()));
			return session;
		}
		if (child.Elapsed() > 60s) {
			session.error = "No render window appeared within 60s";
			return session;
		}
	}
	session.timeToWindow = child.Elapsed();

	const auto dwellEnd = std::chrono::steady_clock::now() + dwell;
	while (std::chrono::steady_clock::now() < dwellEnd) {
		if (auto code = child.WaitFor(250ms)) {
			session.error = Lima::Format("render exited with {} while its window should have stayed open. Output:\n{}", DescribeExitCode(*code), Tail(child.Output()));
			return session;
		}
		DWORD_PTR ignored = 0;
		if (!SendMessageTimeoutW(window, WM_NULL, 0, 0, SMTO_ABORTIFHUNG, 3000, &ignored)) {
			session.error = "The render window stopped responding";
			return session;
		}
	}

	if (screenshot)
		session.distinctColors = CaptureWindow(window, *screenshot);

	PostMessageW(window, WM_CLOSE, 0, 0);
	const auto code = child.WaitFor(15s);
	if (!code)
		session.error = "render did not exit within 15s of its window being closed";
	else if (*code != 0)
		session.error = Lima::Format("render exited with {} after its window was closed. Output:\n{}", DescribeExitCode(*code), Tail(child.Output()));
	return session;
}

// Lets the user see the output of programs without a display of their own. Problems here are only warnings;
// the render test owns render correctness
void Preview(const fs::path& workDir, const std::string& gro, const std::string& top) {
	if (!config.preview) return;
	std::string error;
	try { error = RunRenderAndClose(workDir, { "render", "-f", gro, "-t", top }, config.previewDwell, std::nullopt).error; }
	catch (const TestFailure& failure) { error = failure.message; }
	if (!error.empty()) {
		setConsoleTextColorYellow();
		std::cout << "Preview of " << gro << " failed: " << error << "\n";
		setConsoleTextColorDefault();
	}
}

// ------------------------------------------------ Files ------------------------------------------------ //

struct Vec3 {
	double x = 0, y = 0, z = 0;
	Vec3 operator+(const Vec3& o) const { return { x + o.x, y + o.y, z + o.z }; }
	Vec3 operator-(const Vec3& o) const { return { x - o.x, y - o.y, z - o.z }; }
	Vec3 operator/(double s) const { return { x / s, y / s, z / s }; }
	double Len() const { return std::sqrt(x * x + y * y + z * z); }
	bool Finite() const { return std::isfinite(x) && std::isfinite(y) && std::isfinite(z); }
};

// Minimum image distance in a rectangular periodic box
double PbcDistance(Vec3 a, Vec3 b, Vec3 box) {
	Vec3 d = a - b;
	d.x -= box.x * std::round(d.x / box.x);
	d.y -= box.y * std::round(d.y / box.y);
	d.z -= box.z * std::round(d.z / box.z);
	return d.Len();
}

struct GroAtom {
	int resnr;
	std::string resname;
	std::string name;
	Vec3 position;
};

struct Gro {
	std::string title;
	std::vector<GroAtom> atoms;
	Vec3 box;

	Vec3 Center(size_t begin, size_t end) const {
		Vec3 sum;
		for (size_t i = begin; i < end; i++) sum = sum + atoms[i].position;
		return sum / static_cast<double>(end - begin);
	}

	// Residues as consecutive runs of the same (resnr, resname); resnr wraps at 100000 in .gro files
	size_t ResidueCount() const {
		size_t count = 0;
		for (size_t i = 0; i < atoms.size(); i++)
			if (i == 0 || atoms[i].resnr != atoms[i - 1].resnr || atoms[i].resname != atoms[i - 1].resname) count++;
		return count;
	}
};

std::string Trim(const std::string& s) {
	const auto begin = s.find_first_not_of(" \t\r");
	if (begin == std::string::npos) return "";
	return s.substr(begin, s.find_last_not_of(" \t\r") - begin + 1);
}

// Standalone reader for the standard GROMACS .gro format, so verification does not depend on LIMA's own parser
Gro ReadGro(const fs::path& path) {
	ASSERT(fs::exists(path), Lima::Format("Expected output {} does not exist", path.string()));
	std::ifstream file(path);
	Gro gro;
	std::string line;
	std::getline(file, gro.title);
	ASSERT(std::getline(file, line), Lima::Format("{} has no atom count line", path.string()));
	size_t nAtoms = 0;
	try { nAtoms = std::stoul(Trim(line)); }
	catch (...) { throw TestFailure{ Lima::Format("{} has an invalid atom count line '{}'", path.string(), line) }; }

	gro.atoms.reserve(nAtoms);
	for (size_t i = 0; i < nAtoms; i++) {
		ASSERT(std::getline(file, line), Lima::Format("{} declares {} atoms but has only {}", path.string(), nAtoms, i));
		ASSERT(line.size() > 20, Lima::Format("{} line {} is too short: '{}'", path.string(), i + 3, line));
		GroAtom atom;
		std::istringstream positions(line.substr(20));
		ASSERT(positions >> atom.position.x >> atom.position.y >> atom.position.z,
			Lima::Format("{} line {} has no valid position: '{}'", path.string(), i + 3, line));
		try { atom.resnr = std::stoi(line.substr(0, 5)); }
		catch (...) { throw TestFailure{ Lima::Format("{} line {} has an invalid residue number", path.string(), i + 3) }; }
		atom.resname = Trim(line.substr(5, 5));
		atom.name = Trim(line.substr(10, 5));
		gro.atoms.push_back(std::move(atom));
	}
	ASSERT(std::getline(file, line), Lima::Format("{} has no box line", path.string()));
	std::istringstream box(line);
	ASSERT(box >> gro.box.x >> gro.box.y >> gro.box.z, Lima::Format("{} has an invalid box line '{}'", path.string(), line));
	return gro;
}

void RequireFinite(const Gro& gro, const std::string& name) {
	for (size_t i = 0; i < gro.atoms.size(); i++)
		ASSERT(gro.atoms[i].position.Finite(), Lima::Format("{}: atom {} has a non-finite position", name, i));
}

void RequireBox(const Gro& gro, double size, const std::string& name) {
	ASSERT(std::abs(gro.box.x - size) < 1e-3 && std::abs(gro.box.y - size) < 1e-3 && std::abs(gro.box.z - size) < 1e-3,
		Lima::Format("{}: expected a {} nm cubic box, got {} {} {}", name, size, gro.box.x, gro.box.y, gro.box.z));
}

void RequireSameAtoms(const Gro& a, const Gro& b, const std::string& name) {
	ASSERT(a.atoms.size() == b.atoms.size(), Lima::Format("{}: expected {} atoms, got {}", name, a.atoms.size(), b.atoms.size()));
	for (size_t i = 0; i < a.atoms.size(); i++)
		ASSERT(a.atoms[i].name == b.atoms[i].name && a.atoms[i].resname == b.atoms[i].resname,
			Lima::Format("{}: atom {} changed from {}/{} to {}/{}", name, i, a.atoms[i].resname, a.atoms[i].name, b.atoms[i].resname, b.atoms[i].name));
}

double MaxDisplacement(const Gro& a, const Gro& b) {
	double maxDisplacement = 0;
	for (size_t i = 0; i < a.atoms.size(); i++)
		maxDisplacement = std::max(maxDisplacement, PbcDistance(a.atoms[i].position, b.atoms[i].position, a.box));
	return maxDisplacement;
}

// Molecule name -> number of instances, from the [ molecules ] section
std::map<std::string, int> MoleculeCounts(const fs::path& top) {
	TopologyFile topology{ top };
	std::map<std::string, int> counts;
	for (const auto& entry : topology.GetSystem().molecules) counts[entry.name] += entry.count;
	return counts;
}

// The number of atoms the topology describes. Must always equal the atom count of the matching .gro file
size_t TopologyAtomCount(const fs::path& top) {
	ASSERT(fs::exists(top), Lima::Format("Expected output {} does not exist", top.string()));
	TopologyFile topology{ top };
	ASSERT(topology.HasSystem(), Lima::Format("{} has no [ system ]", top.string()));
	size_t count = 0;
	for (const auto& entry : topology.GetSystem().molecules) {
		ASSERT(entry.moleculetype != nullptr, Lima::Format("{}: molecule {} has no definition", top.string(), entry.name));
		count += entry.moleculetype->atoms.size() * entry.count;
	}
	return count;
}

void RequireTopologyMatches(const fs::path& top, const Gro& gro) {
	const size_t topologyAtoms = TopologyAtomCount(top);
	ASSERT(topologyAtoms == gro.atoms.size(), Lima::Format("{} describes {} atoms, but the coordinates contain {}",
		top.filename().string(), topologyAtoms, gro.atoms.size()));
}

struct Trr {
	int nAtoms = 0;
	int nFrames = 0;
	std::vector<Vec3> lastFrame;
};

// Reads a GROMACS .trr trajectory with the standard xdrfile library
Trr ReadTrr(const fs::path& path) {
	ASSERT(fs::exists(path), Lima::Format("Expected output {} does not exist", path.string()));
	std::string pathString = path.string();
	Trr trr;
	ASSERT(read_trr_natoms(pathString.data(), &trr.nAtoms) == exdrOK && trr.nAtoms > 0, Lima::Format("{} is not a readable .trr file", path.filename().string()));

	XDRFILE* file = xdrfile_open(pathString.c_str(), "r");
	ASSERT(file != nullptr, Lima::Format("Could not open {}", path.string()));
	std::vector<float> x(3 * static_cast<size_t>(trr.nAtoms));
	matrix box;
	int step = 0;
	float time = 0, lambda = 0;
	while (read_trr(file, trr.nAtoms, &step, &time, &lambda, box, reinterpret_cast<rvec*>(x.data()), nullptr, nullptr) == exdrOK)
		trr.nFrames++;
	xdrfile_close(file);

	for (size_t i = 0; i < x.size(); i += 3)
		trr.lastFrame.push_back({ x[i], x[i + 1], x[i + 2] });
	return trr;
}

std::string ReadBytes(const fs::path& path) {
	std::ifstream file(path, std::ios::binary);
	return { std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>() };
}

void CopyFixtures(const fs::path& workDir, std::initializer_list<const char*> files) {
	for (const char* file : files)
		fs::copy_file(config.fixtures / file, workDir / file, fs::copy_options::overwrite_existing);
}

void Banner(const std::string& text) {
	std::cout << "\n\033[36m=== " << text << " ===\033[0m\n" << std::flush;
}

// ------------------------------------------------ Tests ------------------------------------------------ //

TestRoutine TestDispatcher(fs::path workDir) {
	Banner("lima dispatcher");
	RunOptions quiet{ .timeout = 60s, .captureStdout = true, .echo = false };

	auto general = Run(workDir, {}, quiet);
	RequireSuccess(general);
	for (const auto& command : Cli::Commands)
		ASSERT(general.output.find(command.name) != std::string::npos, Lima::Format("'lima' does not list the command '{}'", command.name));

	auto version = Run(workDir, { "--version" }, quiet);
	RequireSuccess(version);
	ASSERT(!Trim(version.output).empty(), "'lima --version' printed nothing");

	for (const auto& command : Cli::Commands) {
		for (const auto& args : { std::vector<std::string>{ "help", std::string(command.name) }, std::vector<std::string>{ std::string(command.name), "--help" } }) {
			auto help = Run(workDir, args, quiet);
			RequireSuccess(help);
			ASSERT(!Trim(help.output).empty(), Lima::Format("'{}' printed nothing", help.commandline));
			ASSERT(!help.windowSeen, Lima::Format("'{}' opened a window", help.commandline));
		}
	}

	RequireExitCode(Run(workDir, { "notacommand" }, quiet), 2);
	RequireExitCode(Run(workDir, { "makebox" }, quiet), 2);	// Missing required --box-size
	RequireExitCode(Run(workDir, { "makebox", "--box-size", "5", "--notanoption" }, quiet), 2);
	ASSERT(fs::is_empty(workDir), "A command that was rejected for invalid arguments still created files");

	co_return LimaUnittestResult{ true, Lima::Format("{} commands respond to help, bad input exits 2", Cli::Commands.size()), false };
}

TestRoutine TestMakeSimParams(fs::path workDir) {
	Banner("lima makesimparams");
	RequireSuccess(Run(workDir, { "makesimparams" }));

	const fs::path path = workDir / "sim_params.txt";
	ASSERT(fs::exists(path), "sim_params.txt was not created");
	std::ifstream file(path);
	std::set<std::string> keys;
	std::string line;
	while (std::getline(file, line)) {
		line = Trim(line.substr(0, line.find('#')));
		if (line.empty() || line.starts_with("//")) continue;
		const auto equals = line.find('=');
		ASSERT(equals != std::string::npos, Lima::Format("sim_params.txt has a line that is not 'key = value': '{}'", line));
		keys.insert(Trim(line.substr(0, equals)));
	}
	ASSERT(keys.contains("n_steps") && keys.contains("dt"), "sim_params.txt does not define n_steps and dt");

	co_return LimaUnittestResult{ true, Lima::Format("{} parameters", keys.size()), false };
}

TestRoutine TestMakeBox(fs::path workDir) {
	Banner("lima makebox");
	RequireSuccess(Run(workDir, { "makebox", "--box-size", "7", "--name", "box" }));

	const Gro gro = ReadGro(workDir / "box.gro");
	ASSERT(gro.atoms.empty(), Lima::Format("box.gro should be empty, but has {} atoms", gro.atoms.size()));
	RequireBox(gro, 7, "box.gro");
	ASSERT(fs::exists(workDir / "box.top"), "box.top was not created");

	co_return LimaUnittestResult{ true, "Empty 7 nm box", false };
}

TestRoutine TestEditConf(fs::path workDir) {
	Banner("lima editconf");
	CopyFixtures(workDir, { "met.gro", "met.top" });
	const std::string inputBytes = ReadBytes(workDir / "met.gro");

	RequireSuccess(Run(workDir, { "editconf", "-c", "met.gro", "-t", "met.top", "--conf-out", "out.gro",
		"--set-center", "1.5", "1.5", "1.5", "--rotate", "0", "0", "1.5708" }));

	ASSERT(ReadBytes(workDir / "met.gro") == inputBytes, "editconf modified its input although --conf-out was given");
	const Gro in = ReadGro(workDir / "met.gro");
	const Gro out = ReadGro(workDir / "out.gro");
	RequireSameAtoms(in, out, "out.gro");
	RequireFinite(out, "out.gro");

	const Vec3 center = out.Center(0, out.atoms.size());
	// Loose, since 'center' may mean the mean position or the bounding box center, and rotating shifts the latter
	ASSERT((center - Vec3{ 1.5, 1.5, 1.5 }).Len() < 0.1, Lima::Format("Center is ({:.3f} {:.3f} {:.3f}), expected (1.5 1.5 1.5)", center.x, center.y, center.z));

	// A rigid transformation preserves every intramolecular distance
	double maxDistanceError = 0;
	for (size_t i = 0; i < in.atoms.size(); i++)
		for (size_t j = i + 1; j < in.atoms.size(); j++)
			maxDistanceError = std::max(maxDistanceError, std::abs(
				(in.atoms[i].position - in.atoms[j].position).Len() - (out.atoms[i].position - out.atoms[j].position).Len()));
	ASSERT(maxDistanceError < 0.005, Lima::Format("The molecule was deformed: an interatomic distance changed by {:.4f} nm", maxDistanceError));

	// Compare after moving both to the same center, so only the rotation remains
	const Vec3 inCenter = in.Center(0, in.atoms.size());
	double maxRotationDisplacement = 0;
	for (size_t i = 0; i < in.atoms.size(); i++)
		maxRotationDisplacement = std::max(maxRotationDisplacement, ((out.atoms[i].position - center) - (in.atoms[i].position - inCenter)).Len());
	ASSERT(maxRotationDisplacement > 0.05, "The molecule was not rotated");

	Preview(workDir, "out.gro", "met.top");
	co_return LimaUnittestResult{ true, Lima::Format("Rigid, centered, rotated (distance error {:.4f} nm)", maxDistanceError), false };
}

TestRoutine TestToGmx(fs::path workDir) {
	Banner("lima togmx");
	CopyFixtures(workDir, { "6lzm.pdb" });
	RequireSuccess(Run(workDir, { "togmx", "-f", "6lzm.pdb", "--name", "lzm" }, { .timeout = 300s }));

	// Distinct protein residues and atoms in the input
	std::set<std::string> pdbResidues;
	std::set<std::string> pdbAtoms;
	{
		std::ifstream pdb(workDir / "6lzm.pdb");
		std::string line;
		while (std::getline(pdb, line)) {
			if (!line.starts_with("ATOM") || line.size() < 27) continue;
			pdbResidues.insert(line.substr(21, 6));	// Chain, residue number and insertion code
			pdbAtoms.insert(line.substr(12, 4) + line.substr(21, 6));
		}
	}

	const Gro gro = ReadGro(workDir / "lzm.gro");
	RequireFinite(gro, "lzm.gro");
	ASSERT(gro.ResidueCount() == pdbResidues.size(), Lima::Format("lzm.gro has {} residues, the input has {}", gro.ResidueCount(), pdbResidues.size()));
	ASSERT(gro.atoms.size() > pdbAtoms.size() && gro.atoms.size() < 3 * pdbAtoms.size(),
		Lima::Format("lzm.gro has {} atoms; with hydrogens added, expected between {} and {}", gro.atoms.size(), pdbAtoms.size(), 3 * pdbAtoms.size()));
	RequireTopologyMatches(workDir / "lzm.top", gro);

	bool hasRestraints = false;
	for (const auto& entry : fs::directory_iterator(workDir))
		hasRestraints |= entry.path().extension() == ".itp" && entry.path().filename().string().find("posre") != std::string::npos;
	ASSERT(hasRestraints, "No position restraint file was written");

	Preview(workDir, "lzm.gro", "lzm.top");
	co_return LimaUnittestResult{ true, Lima::Format("{} residues, {} atoms", gro.ResidueCount(), gro.atoms.size()), false };
}

TestRoutine TestSolvate(fs::path workDir) {
	Banner("lima solvate");
	CopyFixtures(workDir, { "met_box4.gro", "met_box4.top" });
	const std::string inputBytes = ReadBytes(workDir / "met_box4.gro");
	const auto inputMolecules = MoleculeCounts(workDir / "met_box4.top");

	RequireSuccess(Run(workDir, { "solvate", "-c", "met_box4.gro", "-t", "met_box4.top" }));

	ASSERT(ReadBytes(workDir / "met_box4.gro") == inputBytes, "solvate modified its input coordinates");
	const Gro in = ReadGro(workDir / "met_box4.gro");
	const Gro out = ReadGro(workDir / "met_box4_solvated.gro");
	RequireFinite(out, "met_box4_solvated.gro");
	RequireBox(out, 4, "met_box4_solvated.gro");
	RequireTopologyMatches(workDir / "met_box4_solvated.top", out);

	for (size_t i = 0; i < in.atoms.size(); i++)
		ASSERT(out.atoms[i].name == in.atoms[i].name && (out.atoms[i].position - in.atoms[i].position).Len() < 0.002,
			Lima::Format("Solute atom {} was changed by solvate", i));
	for (size_t i = in.atoms.size(); i < out.atoms.size(); i++) {
		const Vec3 p = out.atoms[i].position;
		ASSERT(p.x > -0.2 && p.y > -0.2 && p.z > -0.2 && p.x < out.box.x + 0.2 && p.y < out.box.y + 0.2 && p.z < out.box.z + 0.2,
			Lima::Format("Water atom {} is outside the box", i));
	}

	int solventMolecules = 0;
	for (const auto& [name, count] : MoleculeCounts(workDir / "met_box4_solvated.top"))
		solventMolecules += count - (inputMolecules.contains(name) ? inputMolecules.at(name) : 0);
	ASSERT(solventMolecules > 0, "No water molecules were added to the topology");

	const double expected = SimulationBuilder::defaultSolventsPerNm3 * out.box.x * out.box.y * out.box.z;
	ASSERT(solventMolecules > 0.5 * expected && solventMolecules < 1.2 * expected,
		Lima::Format("Added {} waters, expected about {:.0f} at the default density of {}/nm^3",
			solventMolecules, expected, SimulationBuilder::defaultSolventsPerNm3));

	Preview(workDir, "met_box4_solvated.gro", "met_box4_solvated.top");
	co_return LimaUnittestResult{ true, Lima::Format("{} waters added", solventMolecules), false };
}

TestRoutine TestInsertMolecule(fs::path workDir) {
	Banner("lima insertmolecule");
	CopyFixtures(workDir, { "met.gro", "met.top", "box5.gro", "box5.top" });
	const size_t moleculeAtoms = ReadGro(workDir / "met.gro").atoms.size();

	RequireSuccess(Run(workDir, { "insertmolecule", "--conf-source", "met.gro", "--top-source", "met.top",
		"--conf-target", "box5.gro", "--top-target", "box5.top", "--position", "2", "2", "2" }));

	const Gro out = ReadGro(workDir / "box5.gro");
	ASSERT(out.atoms.size() == moleculeAtoms, Lima::Format("The target box has {} atoms after inserting one {}-atom molecule", out.atoms.size(), moleculeAtoms));
	RequireFinite(out, "box5.gro");
	RequireBox(out, 5, "box5.gro");
	RequireTopologyMatches(workDir / "box5.top", out);
	const Vec3 center = out.Center(0, out.atoms.size());
	ASSERT((center - Vec3{ 2, 2, 2 }).Len() < 0.1, Lima::Format("The molecule was inserted at ({:.2f} {:.2f} {:.2f}), expected (2 2 2)", center.x, center.y, center.z));

	Preview(workDir, "box5.gro", "box5.top");
	co_return LimaUnittestResult{ true, "Inserted at the requested position", false };
}

TestRoutine TestInsertMolecules(fs::path workDir) {
	Banner("lima insertmolecules");
	CopyFixtures(workDir, { "met.gro", "met.top", "met_box6.gro", "met_box6.top" });
	const size_t moleculeAtoms = ReadGro(workDir / "met.gro").atoms.size();
	const size_t targetAtoms = ReadGro(workDir / "met_box6.gro").atoms.size();
	constexpr int nInsertions = 10;

	auto result = Run(workDir, { "insertmolecules", "--conf-source", "met.gro", "--top-source", "met.top",
		"--conf-target", "met_box6.gro", "--top-target", "met_box6.top", "-n", std::to_string(nInsertions), "--rotate-randomly", "-d" },
		{ .timeout = 600s });
	RequireSuccess(result);
	RequireWindowSeen(result);

	const Gro out = ReadGro(workDir / "met_box6.gro");
	const size_t expectedAtoms = targetAtoms + nInsertions * moleculeAtoms;
	ASSERT(out.atoms.size() == expectedAtoms, Lima::Format("Expected {} atoms after inserting {} molecules, got {}", expectedAtoms, nInsertions, out.atoms.size()));
	RequireFinite(out, "met_box6.gro");
	RequireBox(out, 6, "met_box6.gro");
	RequireTopologyMatches(workDir / "met_box6.top", out);

	// Every molecule is a consecutive block of atoms. No two molecules should be on top of each other
	double minDistance = DBL_MAX;
	for (size_t i = 0; i < out.atoms.size(); i++)
		for (size_t j = (i / moleculeAtoms + 1) * moleculeAtoms; j < out.atoms.size(); j++)
			minDistance = std::min(minDistance, PbcDistance(out.atoms[i].position, out.atoms[j].position, out.box));
	ASSERT(minDistance > 0.08, Lima::Format("Two molecules overlap: atoms only {:.3f} nm apart", minDistance));

	co_return LimaUnittestResult{ true, Lima::Format("{} molecules, closest contact {:.2f} nm", nInsertions, minDistance), false };
}

TestRoutine TestEnergyMinimization(fs::path workDir) {
	Banner("lima em");
	// Two waters in this fixture overlap; minimization must push them apart
	CopyFixtures(workDir, { "metsol_clash.gro", "metsol.top", "Protein_chain_A.itp", "SOL.itp" });
	constexpr size_t clashA = 19, clashB = 22;	// Oxygens of the overlapping waters
	const std::string inputBytes = ReadBytes(workDir / "metsol_clash.gro");

	auto result = Run(workDir, { "em", "-c", "metsol_clash.gro", "-t", "metsol.top", "--conf-out", "em.gro", "-d" }, { .timeout = 600s });
	RequireSuccess(result);
	RequireWindowSeen(result);

	ASSERT(ReadBytes(workDir / "metsol_clash.gro") == inputBytes, "em modified its input although --conf-out was given");
	const Gro in = ReadGro(workDir / "metsol_clash.gro");
	const Gro out = ReadGro(workDir / "em.gro");
	RequireSameAtoms(in, out, "em.gro");
	RequireFinite(out, "em.gro");

	const double before = PbcDistance(in.atoms[clashA].position, in.atoms[clashB].position, in.box);
	const double after = PbcDistance(out.atoms[clashA].position, out.atoms[clashB].position, out.box);
	ASSERT(after > 0.2, Lima::Format("The overlapping waters were not separated: {:.3f} nm -> {:.3f} nm", before, after));

	co_return LimaUnittestResult{ true, Lima::Format("Clash resolved: {:.2f} nm -> {:.2f} nm", before, after), false };
}

TestRoutine TestMdrun(fs::path workDir) {
	Banner("lima mdrun");
	CopyFixtures(workDir, { "metsol.gro", "metsol.top", "Protein_chain_A.itp", "SOL.itp", "md_params.txt" });

	auto result = Run(workDir, { "mdrun", "-c", "metsol.gro", "-t", "metsol.top", "-s", "md_params.txt",
		"--conf-out", "out.gro", "--trajectory", "traj.trr", "--uff", "-d" }, { .timeout = 900s });
	RequireSuccess(result);
	RequireWindowSeen(result);

	const Gro in = ReadGro(workDir / "metsol.gro");
	const Gro out = ReadGro(workDir / "out.gro");
	RequireSameAtoms(in, out, "out.gro");
	RequireFinite(out, "out.gro");
	for (size_t i = 0; i < out.atoms.size(); i++) {
		const Vec3 p = out.atoms[i].position;
		ASSERT(p.x > -1 && p.y > -1 && p.z > -1 && p.x < out.box.x + 1 && p.y < out.box.y + 1 && p.z < out.box.z + 1,
			Lima::Format("Atom {} ended far outside the box", i));
	}
	const double moved = MaxDisplacement(in, out);
	ASSERT(moved > 0.01, "No atom moved during the simulation");

	// md_params.txt runs 20000 steps and logs every 200th
	constexpr int maxFrames = 20000 / 200 + 1;
	const Trr trr = ReadTrr(workDir / "traj.trr");
	ASSERT(trr.nAtoms == static_cast<int>(in.atoms.size()), Lima::Format("The trajectory has {} atoms, the system has {}", trr.nAtoms, in.atoms.size()));
	ASSERT(trr.nFrames >= 2 && trr.nFrames <= maxFrames,
		Lima::Format("The trajectory has {} frames, expected at most {} with the requested logging interval", trr.nFrames, maxFrames));
	// The last logged frame is at most one logging interval (0.4 ps) before the final coordinates
	double meanDistanceToFinal = 0;
	for (size_t i = 0; i < out.atoms.size(); i++) {
		ASSERT(trr.lastFrame[i].Finite(), Lima::Format("Atom {} has a non-finite position in the trajectory", i));
		meanDistanceToFinal += PbcDistance(trr.lastFrame[i], out.atoms[i].position, out.box) / out.atoms.size();
	}
	ASSERT(meanDistanceToFinal < 0.2, Lima::Format("The last trajectory frame does not resemble the final coordinates (mean distance {:.2f} nm)", meanDistanceToFinal));
	const auto trajectoryBytes = fs::file_size(workDir / "traj.trr");
	ASSERT(fs::exists(workDir / "trajectory.uff") && fs::file_size(workDir / "trajectory.uff") > 0, "--uff did not write trajectory.uff");

	co_return LimaUnittestResult{ true, Lima::Format("{} atoms, trajectory {:.1f} MB", out.atoms.size(), trajectoryBytes / 1e6), false };
}

TestRoutine TestBuildMembrane(fs::path workDir) {
	Banner("lima buildmembrane");
	constexpr double boxSize = 8;
	auto result = Run(workDir, { "buildmembrane", "--lipids", "POPC", "70", "cholesterol", "30",
		"--box-size", "8", "--seed", "1", "-d" }, { .timeout = 900s });
	RequireSuccess(result);
	RequireWindowSeen(result);

	const Gro gro = ReadGro(workDir / "membrane.gro");
	RequireFinite(gro, "membrane.gro");
	RequireBox(gro, boxSize, "membrane.gro");
	RequireTopologyMatches(workDir / "membrane.top", gro);

	int popc = 0, total = 0;
	for (const auto& [name, count] : MoleculeCounts(workDir / "membrane.top")) {
		total += count;
		if (name.find("POPC") != std::string::npos || name.find("popc") != std::string::npos) popc += count;
	}
	ASSERT(total > 0, "The membrane contains no lipids");
	const double popcFraction = static_cast<double>(popc) / total;
	ASSERT(popcFraction > 0.6 && popcFraction < 0.8, Lima::Format("{:.0f}% of the lipids are POPC, expected 70%", popcFraction * 100));

	const double areaPerLipid = boxSize * boxSize / (total / 2.0);
	ASSERT(areaPerLipid > 0.3 && areaPerLipid < 1.5, Lima::Format("{} lipids give {:.2f} nm^2 per lipid, which is not a plausible bilayer", total, areaPerLipid));

	const double meanZ = gro.Center(0, gro.atoms.size()).z;
	ASSERT(std::abs(meanZ - boxSize / 2) < 0.3, Lima::Format("The membrane is centered at z={:.2f}, expected the default of {}", meanZ, boxSize / 2));

	co_return LimaUnittestResult{ true, Lima::Format("{} lipids, {:.0f}% POPC, {:.2f} nm^2/lipid", total, popcFraction * 100, areaPerLipid), false };
}

TestRoutine TestRender(fs::path workDir) {
	Banner("lima render");
	CopyFixtures(workDir, { "metsol.gro", "metsol.top", "Protein_chain_A.itp", "SOL.itp" });

	const fs::path screenshot = workDir / "render.bmp";
	const auto session = RunRenderAndClose(workDir, { "render", "-f", "metsol.gro", "-t", "metsol.top", "--highlight", "0", "1", "2" },
		config.renderDwell, screenshot);
	ASSERT(session.error.empty(), session.error);
	// Only proves something was drawn, not what; render.bmp is kept for a human to look at
	ASSERT(session.distinctColors > 16, Lima::Format("The render window looks blank ({} distinct colors), see {}", session.distinctColors, screenshot.string()));

	co_return LimaUnittestResult{ true, Lima::Format("Window after {:.1f}s, closed cleanly, screenshot {}",
		session.timeToWindow.count(), screenshot.filename().string()), false };
}

TestRoutine TestSelfTest(fs::path workDir) {
	Banner("lima selftest");
	auto result = Run(workDir, { "selftest" }, { .timeout = 1800s });
	RequireSuccess(result);
	RequireWindowSeen(result);
	co_return LimaUnittestResult{ true, "", false };
}

struct CliTest {
	std::string command;	// The lima command this test covers, or "dispatcher"
	TestRoutine(*function)(fs::path);
};

// Cheapest first. The dispatcher test runs first since every other test depends on lima starting at all
const std::vector<CliTest> cliTests{
	{ "dispatcher", TestDispatcher },
	{ "makesimparams", TestMakeSimParams },
	{ "makebox", TestMakeBox },
	{ "editconf", TestEditConf },
	{ "togmx", TestToGmx },
	{ "solvate", TestSolvate },
	{ "insertmolecule", TestInsertMolecule },
	{ "insertmolecules", TestInsertMolecules },
	{ "em", TestEnergyMinimization },
	{ "mdrun", TestMdrun },
	{ "buildmembrane", TestBuildMembrane },
	{ "render", TestRender },
	{ "selftest", TestSelfTest },
};

// ------------------------------------------------ Setup ------------------------------------------------ //

fs::path DefaultLimaExecutable() {
	wchar_t buffer[MAX_PATH];
	GetModuleFileNameW(nullptr, buffer, MAX_PATH);
	// limaclitest.exe lives in <build>/code/LIMA_TESTS, lima.exe in <build>/code/LIMA
	return fs::path(buffer).parent_path().parent_path() / "LIMA" / "lima.exe";
}

void PrintUsage() {
	std::cout << "Usage: limaclitest [--lima PATH] [--no-preview] [--dwell SECONDS] [COMMAND...]\n";
}

void ParseArguments(int argc, char** argv) {
	config.lima = DefaultLimaExecutable();
	for (int i = 1; i < argc; i++) {
		const std::string arg = argv[i];
		if (arg == "--lima" && i + 1 < argc) config.lima = argv[++i];
		else if (arg == "--no-preview") config.preview = false;
		else if (arg == "--dwell" && i + 1 < argc) config.renderDwell = std::chrono::seconds(std::stoi(argv[++i]));
		else if (arg == "--help" || arg == "-h") { PrintUsage(); std::exit(0); }
		else if (arg.starts_with("-")) { PrintUsage(); throw std::runtime_error("Unknown option " + arg); }
		else {
			if (std::ranges::none_of(cliTests, [&](const CliTest& test) { return test.command == arg; }))
				throw std::runtime_error("No test for command " + arg);
			config.only.insert(arg);
		}
	}
	config.lima = fs::absolute(config.lima);
}

// Old runs are deleted, the current one is kept for inspection
void PrepareRunDir() {
	const fs::path runsDir = FileUtils::GetLimaDir() / "tests" / "clitests" / "_runs";
	if (fs::exists(runsDir))
		for (const auto& entry : fs::directory_iterator(runsDir)) {
			std::error_code error;
			fs::remove_all(entry.path(), error);
			if (error) std::cout << "Could not delete old run " << entry.path().string() << ": " << error.message() << "\n";
		}
	const auto now = std::chrono::floor<std::chrono::seconds>(std::chrono::system_clock::now());
	config.runDir = runsDir / Lima::Format("{:%Y%m%d-%H%M%S}", std::chrono::zoned_time{ std::chrono::current_zone(), now });
	fs::create_directories(config.runDir);
}

void SetupConsoleAndJob() {
	// lima prints ANSI sequences; make sure the console interprets them
	HANDLE out = GetStdHandle(STD_OUTPUT_HANDLE);
	DWORD mode = 0;
	if (GetConsoleMode(out, &mode))
		SetConsoleMode(out, mode | ENABLE_VIRTUAL_TERMINAL_PROCESSING);

	jobObject = CreateJobObjectW(nullptr, nullptr);
	JOBOBJECT_EXTENDED_LIMIT_INFORMATION limits{};
	limits.BasicLimitInformation.LimitFlags = JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE;
	SetInformationJobObject(jobObject, JobObjectExtendedLimitInformation, &limits, sizeof(limits));
}

} // namespace

int main(int argc, char** argv) {
	try {
		SetupConsoleAndJob();
		ParseArguments(argc, argv);
		if (!fs::exists(config.lima))
			throw std::runtime_error("lima executable not found: " + config.lima.string() + " (use --lima PATH)");
		PrepareRunDir();

		std::cout << "Testing " << config.lima.string() << "\n"
			<< "Output in " << config.runDir.string() << "\n"
			<< "Leave the windows alone, they close by themselves\n";

		LimaUnittestManager testman{ false };

		// Every lima command must be covered, so a new command cannot ship untested
		std::vector<std::string> untested;
		for (const auto& command : Cli::Commands)
			if (std::ranges::none_of(cliTests, [&](const CliTest& test) { return test.command == command.name; }))
				untested.emplace_back(command.name);
		testman.AddTest("Every command has a test", [untested]() -> TestRoutine {
			std::string missing;
			for (const auto& name : untested) missing += " " + name;
			co_return LimaUnittestResult{ untested.empty(), untested.empty() ? "" : "No test for:" + missing, false };
		});

		for (const auto& test : cliTests) {
			if (!config.only.empty() && !config.only.contains(test.command)) continue;
			const fs::path workDir = config.runDir / test.command;
			fs::create_directories(workDir);
			testman.AddTest("lima " + test.command, [&test, workDir] { return test.function(workDir); });
		}

		return testman.Finish() == 0 ? 0 : 1;
	}
	catch (const std::exception& ex) {
		std::cerr << "limaclitest: " << ex.what() << std::endl;
		return 2;
	}
}
