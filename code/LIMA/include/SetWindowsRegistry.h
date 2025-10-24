#ifdef _WIN32

#include <windows.h>
#include <string>
#include <shlobj.h> // for SHChangeNotify

namespace {
    bool RegistryValueEquals(HKEY root, const wchar_t* path, const wchar_t* value) {
        HKEY key;
        if (RegOpenKeyExW(root, path, 0, KEY_READ, &key) != ERROR_SUCCESS)
            return false;

        wchar_t buffer[512];
        DWORD size = sizeof(buffer);
        auto status = RegQueryValueExW(key, nullptr, nullptr, nullptr, reinterpret_cast<LPBYTE>(buffer), &size);
        RegCloseKey(key);
        if (status != ERROR_SUCCESS)
            return false;

        return wcscmp(buffer, value) == 0;
    }
}

void RegisterGrofileAssociation() {
    const wchar_t* progID = L"LIMA.grofile";
    const wchar_t* desc = L"LIMA .gro File";
    wchar_t exePath[MAX_PATH];
    GetModuleFileNameW(nullptr, exePath, MAX_PATH);

    // Build command string
    std::wstring command = L"\"" + std::wstring(exePath) + L"\" \"render\" \"-conf\" \"%1\"";

    // Already correct? -> skip
    if (RegistryValueEquals(HKEY_CURRENT_USER, L"Software\\Classes\\.gro", progID) &&
        RegistryValueEquals(HKEY_CURRENT_USER, L"Software\\Classes\\LIMA.grofile\\shell\\open\\command", command.c_str()))
        return;

    HKEY key;

    // 1. Link extension -> progID
    RegCreateKeyW(HKEY_CURRENT_USER, L"Software\\Classes\\.gro", &key);
    RegSetValueW(key, nullptr, REG_SZ, progID, 0);
    RegCloseKey(key);

    // 2. Define progID
    RegCreateKeyW(HKEY_CURRENT_USER, L"Software\\Classes\\LIMA.grofile", &key);
    RegSetValueW(key, nullptr, REG_SZ, desc, 0);
    RegCloseKey(key);

    // 3. Add open command
    RegCreateKeyW(HKEY_CURRENT_USER, L"Software\\Classes\\LIMA.grofile\\shell\\open\\command", &key);
    RegSetValueW(key, nullptr, REG_SZ, command.c_str(), 0);
    RegCloseKey(key);

    // Add icon
    RegCreateKeyW(HKEY_CURRENT_USER, L"Software\\Classes\\LIMA.grofile\\DefaultIcon", &key);
    // “,0” = first icon in the exe
    std::wstring iconPath = std::wstring(exePath) + L",0";
    RegSetValueW(key, nullptr, REG_SZ, iconPath.c_str(), 0);
    RegCloseKey(key);


    // Also register under Applications for Open With menu icon
    {
        wchar_t exePath[MAX_PATH];
        GetModuleFileNameW(nullptr, exePath, MAX_PATH);

        const wchar_t* appKeyPath = L"Software\\Classes\\Applications\\lima.exe";

        HKEY key;
        if (RegCreateKeyW(HKEY_CURRENT_USER, (std::wstring(appKeyPath) + L"\\shell\\open\\command").c_str(), &key) == ERROR_SUCCESS) {
            std::wstring cmd = L"\"" + std::wstring(exePath) + L"\" \"%1\"";
            RegSetValueW(key, nullptr, REG_SZ, cmd.c_str(), 0);
            RegCloseKey(key);
        }

        if (RegCreateKeyW(HKEY_CURRENT_USER, (std::wstring(appKeyPath) + L"\\DefaultIcon").c_str(), &key) == ERROR_SUCCESS) {
            std::wstring icon = std::wstring(exePath) + L",0";
            RegSetValueW(key, nullptr, REG_SZ, icon.c_str(), 0);
            RegCloseKey(key);
        }
    }



    // 4. Refresh Explorer
    SHChangeNotify(SHCNE_ASSOCCHANGED, SHCNF_IDLIST, nullptr, nullptr);
}
#else
void RegisterGrofileAssociation() {
    // No-op on non-Windows platforms
}
#endif

