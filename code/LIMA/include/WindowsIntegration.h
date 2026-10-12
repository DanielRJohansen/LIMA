#pragma once
// Integrates a lima.exe with Windows: its folder on the user's PATH, and .gro files opening with it.
// Everything is under HKEY_CURRENT_USER, so no admin rights are needed. The first interactive run of a release
// lima.exe offers both (FirstRunSetup), and 'lima setup' does both without asking.

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#ifdef _WIN32

#include <windows.h>
#include <shlobj.h> // for SHChangeNotify
#include <cstdio>
#include <cwctype>
#include <io.h>
#include <optional>
#include <stdexcept>

namespace WindowsIntegration {
namespace detail {
    namespace fs = std::filesystem;

    constexpr const wchar_t* environmentKey = L"Environment";
    constexpr const wchar_t* settingsKey = L"Software\\LIMA";
    constexpr const wchar_t* progID = L"LIMA.grofile";

    inline fs::path ExePath() {
        wchar_t buffer[MAX_PATH];
        const DWORD length = GetModuleFileNameW(nullptr, buffer, MAX_PATH);
        return length > 0 && length < MAX_PATH ? fs::path(buffer) : fs::path{};
    }

    inline std::optional<std::wstring> ReadString(HKEY root, const wchar_t* key, const wchar_t* value, DWORD* type = nullptr) {
        // RRF_NOEXPAND keeps %VARIABLES% in PATH, so writing it back does not freeze them
        constexpr DWORD flags = RRF_RT_REG_SZ | RRF_RT_REG_EXPAND_SZ | RRF_NOEXPAND;
        DWORD size = 0;
        if (RegGetValueW(root, key, value, flags, type, nullptr, &size) != ERROR_SUCCESS)
            return std::nullopt;
        std::wstring result(size / sizeof(wchar_t), L'\0');
        if (RegGetValueW(root, key, value, flags, type, result.data(), &size) != ERROR_SUCCESS)
            return std::nullopt;
        result.resize(wcslen(result.c_str()));
        return result;
    }

    inline bool WriteString(HKEY root, const wchar_t* key, const wchar_t* value, const std::wstring& data, DWORD type = REG_SZ) {
        HKEY handle;
        if (RegCreateKeyExW(root, key, 0, nullptr, 0, KEY_SET_VALUE, nullptr, &handle, nullptr) != ERROR_SUCCESS)
            return false;
        const auto status = RegSetValueExW(handle, value, 0, type, reinterpret_cast<const BYTE*>(data.c_str()),
            static_cast<DWORD>((data.size() + 1) * sizeof(wchar_t)));
        RegCloseKey(handle);
        return status == ERROR_SUCCESS;
    }

    // ---------------------------------------------------------------- PATH

    inline std::vector<std::wstring> SplitPath(const std::wstring& path) {
        std::vector<std::wstring> entries;
        size_t start = 0;
        while (start <= path.size()) {
            const size_t end = (std::min)(path.find(L';', start), path.size());
            if (end > start)
                entries.push_back(path.substr(start, end - start));
            start = end + 1;
        }
        return entries;
    }

    // An entry's directory, for comparing: variables expanded, backslashes, lowercase, no trailing separator
    inline std::wstring Normalized(const std::wstring& entry) {
        wchar_t expanded[32768];
        const DWORD length = ExpandEnvironmentStringsW(entry.c_str(), expanded, static_cast<DWORD>(std::size(expanded)));
        std::wstring result = length > 0 && length <= std::size(expanded) ? std::wstring(expanded) : entry;
        for (auto& c : result)
            c = c == L'/' ? L'\\' : static_cast<wchar_t>(std::towlower(c));
        while (result.size() > 3 && result.back() == L'\\')
            result.pop_back();
        return result;
    }

    inline bool Contains(const std::vector<std::wstring>& entries, const fs::path& dir) {
        const std::wstring target = Normalized(dir.wstring());
        for (const auto& entry : entries)
            if (Normalized(entry) == target)
                return true;
        return false;
    }

    // A release: resources next to the exe, and not inside the repository. Dev builds also get resources installed
    // next to them by CMake, but are the developer's own business
    inline bool IsReleaseInstall(const fs::path& dir) {
        std::error_code error;
        if (!fs::exists(dir / "resources", error))
            return false;
        for (fs::path path = dir; !path.empty() && path != path.parent_path(); path = path.parent_path())
            if (fs::exists(path / "lima_main_dir.txt", error))
                return false;
        return true;
    }

    // A folder holding another LIMA release, typically an older one, or one that was deleted after unzipping
    // a newer release
    inline bool IsOtherLimaEntry(const std::wstring& entry, const fs::path& dir) {
        if (Normalized(entry) == Normalized(dir.wstring()))
            return false;
        wchar_t expanded[32768];
        const DWORD length = ExpandEnvironmentStringsW(entry.c_str(), expanded, static_cast<DWORD>(std::size(expanded)));
        const fs::path path = length > 0 && length <= std::size(expanded) ? fs::path(expanded) : fs::path(entry);
        std::error_code error;
        if (fs::exists(path / "lima.exe", error) && IsReleaseInstall(path))
            return true;
        const std::wstring folder = Normalized(path.filename().wstring());
        return !fs::exists(path, error) && folder.starts_with(L"lima-") && folder.ends_with(L"-windows-x64");
    }

    inline void SetUserPath(const std::vector<std::wstring>& entries, DWORD type) {
        std::wstring joined;
        for (const auto& entry : entries)
            joined += (joined.empty() ? L"" : L";") + entry;
        if (!WriteString(HKEY_CURRENT_USER, environmentKey, L"Path", joined, type))
            throw std::runtime_error("Could not write the user PATH in the registry");
        // Tells Explorer, and the terminals it starts from now on, about the new PATH
        SendMessageTimeoutW(HWND_BROADCAST, WM_SETTINGCHANGE, 0, reinterpret_cast<LPARAM>(environmentKey),
            SMTO_ABORTIFHUNG, 5000, nullptr);
    }

    inline bool OnPath(const fs::path& dir) {
        if (const auto userPath = ReadString(HKEY_CURRENT_USER, environmentKey, L"Path"); userPath && Contains(SplitPath(*userPath), dir))
            return true;
        // Also covers the system PATH, and folders added by hand to this terminal's PATH
        const wchar_t* processPath = _wgetenv(L"PATH");
        return processPath && Contains(SplitPath(processPath), dir);
    }

    inline std::vector<std::wstring> OtherLimaEntries(const fs::path& dir) {
        std::vector<std::wstring> others;
        if (const auto userPath = ReadString(HKEY_CURRENT_USER, environmentKey, L"Path"))
            for (const auto& entry : SplitPath(*userPath))
                if (IsOtherLimaEntry(entry, dir))
                    others.push_back(entry);
        return others;
    }

    // Adds dir to the user PATH, and removes other copies of LIMA from it
    inline void AddToPath(const fs::path& dir) {
        DWORD type = REG_EXPAND_SZ;
        const auto userPath = ReadString(HKEY_CURRENT_USER, environmentKey, L"Path", &type);
        std::vector<std::wstring> entries;
        for (const auto& entry : SplitPath(userPath.value_or(L"")))
            if (!IsOtherLimaEntry(entry, dir))
                entries.push_back(entry);
        if (!Contains(entries, dir))
            entries.push_back(dir.wstring());
        SetUserPath(entries, type);
    }

    inline bool RemoveFromPath(const fs::path& dir) {
        DWORD type = REG_EXPAND_SZ;
        const auto userPath = ReadString(HKEY_CURRENT_USER, environmentKey, L"Path", &type);
        if (!userPath)
            return false;
        const std::wstring target = Normalized(dir.wstring());
        std::vector<std::wstring> entries;
        for (const auto& entry : SplitPath(*userPath))
            if (Normalized(entry) != target)
                entries.push_back(entry);
        if (entries.size() == SplitPath(*userPath).size())
            return false;
        SetUserPath(entries, type);
        return true;
    }

    // ---------------------------------------------------------------- .gro association

    inline std::wstring OpenCommand(const fs::path& exe) {
        return L"\"" + exe.wstring() + L"\" \"render\" \"-conf\" \"%1\"";
    }

    inline bool GrofileAssociationIsCurrent(const fs::path& exe) {
        return ReadString(HKEY_CURRENT_USER, L"Software\\Classes\\.gro", nullptr) == progID
            && ReadString(HKEY_CURRENT_USER, L"Software\\Classes\\LIMA.grofile\\shell\\open\\command", nullptr) == OpenCommand(exe);
    }

    // ---------------------------------------------------------------- First run

    inline bool IsInteractive() {
        return _isatty(_fileno(stdin)) && _isatty(_fileno(stdout));
    }

    inline bool AskYesNo(const std::string& question) {
        std::cout << question << " [Y/n] " << std::flush;
        std::string answer;
        if (!std::getline(std::cin, answer))
            return false;
        const auto first = answer.find_first_not_of(" \t");
        return first == std::string::npos || answer[first] == 'y' || answer[first] == 'Y';
    }
}

inline void RegisterGrofileAssociation() {
    const std::filesystem::path exe = detail::ExePath();
    if (detail::GrofileAssociationIsCurrent(exe))
        return;
    const std::wstring exePath = exe.wstring();
    const std::wstring icon = exePath + L",0";	// The first icon in the exe

    detail::WriteString(HKEY_CURRENT_USER, L"Software\\Classes\\.gro", nullptr, detail::progID);
    detail::WriteString(HKEY_CURRENT_USER, L"Software\\Classes\\LIMA.grofile", nullptr, L"LIMA .gro File");
    detail::WriteString(HKEY_CURRENT_USER, L"Software\\Classes\\LIMA.grofile\\shell\\open\\command", nullptr, detail::OpenCommand(exe));
    detail::WriteString(HKEY_CURRENT_USER, L"Software\\Classes\\LIMA.grofile\\DefaultIcon", nullptr, icon);

    // Also register under Applications, for the icon in the Open With menu
    detail::WriteString(HKEY_CURRENT_USER, L"Software\\Classes\\Applications\\lima.exe\\shell\\open\\command", nullptr, L"\"" + exePath + L"\" \"%1\"");
    detail::WriteString(HKEY_CURRENT_USER, L"Software\\Classes\\Applications\\lima.exe\\DefaultIcon", nullptr, icon);

    SHChangeNotify(SHCNE_ASSOCCHANGED, SHCNF_IDLIST, nullptr, nullptr);	// Refresh Explorer
}

inline void UnregisterGrofileAssociation() {
    // Other programs may have registered themselves under .gro too, so only LIMA's default is removed
    if (detail::ReadString(HKEY_CURRENT_USER, L"Software\\Classes\\.gro", nullptr) == detail::progID)
        RegDeleteKeyValueW(HKEY_CURRENT_USER, L"Software\\Classes\\.gro", nullptr);
    RegDeleteTreeW(HKEY_CURRENT_USER, L"Software\\Classes\\LIMA.grofile");
    RegDeleteTreeW(HKEY_CURRENT_USER, L"Software\\Classes\\Applications\\lima.exe");
    SHChangeNotify(SHCNE_ASSOCCHANGED, SHCNF_IDLIST, nullptr, nullptr);
}

// 'lima setup': puts this lima.exe on the PATH and opens .gro files with it, or with undo removes both
inline void Setup(std::ostream& out, bool undo) {
    const std::filesystem::path exe = detail::ExePath();
    const std::filesystem::path dir = exe.parent_path();
    if (undo) {
        out << (detail::RemoveFromPath(dir) ? "Removed " + dir.string() + " from your PATH\n" : dir.string() + " was not on your PATH\n");
        UnregisterGrofileAssociation();
        out << ".gro files no longer open with LIMA\n";
        return;
    }
    for (const auto& other : detail::OtherLimaEntries(dir))
        out << "Removed " << std::filesystem::path(other).string() << " from your PATH\n";
    detail::AddToPath(dir);
    out << "Added " << dir.string() << " to your PATH. Terminals opened from now on will find 'lima'\n";
    RegisterGrofileAssociation();
    out << ".gro files now open with LIMA\n";
}

// On the first interactive run of a release lima.exe, offers to put it on the PATH and to open .gro files with it.
// Asked once per lima.exe location, so unzipping a new release asks again, and says nothing if both are already
// in place. Runs from scripts, pipes and tests never ask, as stdin or stdout is not a console there.
inline void FirstRunSetup() {
    const std::filesystem::path exe = detail::ExePath();
    const std::filesystem::path dir = exe.parent_path();
    if (exe.empty() || !detail::IsInteractive() || !detail::IsReleaseInstall(dir))
        return;
    const bool offerPath = !detail::OnPath(dir);
    const bool offerAssociation = !detail::GrofileAssociationIsCurrent(exe);
    if (!offerPath && !offerAssociation)
        return;
    if (detail::ReadString(HKEY_CURRENT_USER, detail::settingsKey, L"SetupAskedFor") == exe.wstring())
        return;

    try {
        std::cout << "LIMA in " << dir.string() << " is not set up for this user yet.\n";
        if (offerPath) {
            const auto others = detail::OtherLimaEntries(dir);
            for (const auto& other : others)
                std::cout << "  Your PATH has another LIMA in " << std::filesystem::path(other).string() << ", which this replaces\n";
            if (detail::AskYesNo("Add this folder to your PATH, so 'lima' works in any terminal?")) {
                detail::AddToPath(dir);
                std::cout << "  Done. Terminals opened from now on will find 'lima'\n";
            }
        }
        if (offerAssociation && detail::AskYesNo("Open .gro files with LIMA when double-clicked?"))
            RegisterGrofileAssociation();
        std::cout << "Run 'lima setup' or 'lima setup --undo' to change this later.\n\n";
    }
    catch (const std::exception& error) {
        std::cout << "LIMA setup failed: " << error.what() << "\n\n";
    }
    detail::WriteString(HKEY_CURRENT_USER, detail::settingsKey, L"SetupAskedFor", exe.wstring());
}

// True when lima was started by double-clicking it in Explorer, so its console window closes as soon as it exits
inline bool OwnsConsoleWindow() {
    DWORD processes[2];
    return detail::IsInteractive() && GetConsoleProcessList(processes, 2) == 1;
}

} // namespace WindowsIntegration

#else

namespace WindowsIntegration {
    inline void RegisterGrofileAssociation() {}
    inline void FirstRunSetup() {}
    inline bool OwnsConsoleWindow() { return false; }
}

#endif
