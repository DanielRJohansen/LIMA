@echo off
rem Builds, tests and packages a LIMA release for Windows and Linux, and uploads it as a draft GitHub release.
rem
rem   distribution\release.bat [--dry-run | --upload-only]
rem
rem It shows the version in the top-level CMakeLists.txt and asks which version to release. Choosing a new one
rem bumps CMakeLists.txt, and commits and pushes that change, before anything is built.
rem --dry-run builds, tests and packages everything into distribution\out, but allows uncommitted changes
rem and skips the upload. It always uses the current version, and the Linux build uses the last commit.
rem --upload-only skips the builds and uploads the files already in distribution\out\<version>, for retrying an
rem upload that failed. It releases the commit those files were built from.
rem
rem One-time setup: see distribution\README.md
setlocal EnableDelayedExpansion

set "VS_DIR=C:\Program Files\Microsoft Visual Studio\2022\Community"
set "WSL_DISTRO=Ubuntu-24.04"
set "GITHUB_REPO=DanielRJohansen/LIMA"
rem RTX 40 (89), H100 (90), B200 (100), RTX 50 (120), plus PTX of the newest. Keep in sync with linux\config.sh
set "CUDA_ARCHITECTURES=89-real;90-real;100-real;120"

set DRY_RUN=0
set UPLOAD_ONLY=0
if "%~1"=="--dry-run" set DRY_RUN=1
if "%~1"=="--upload-only" set UPLOAD_ONLY=1
if not "%~1"=="" if not "%~1"=="--dry-run" if not "%~1"=="--upload-only" (
    echo Usage: release.bat [--dry-run ^| --upload-only]
    exit /b 2
)

pushd "%~dp0.."
set "REPO=%CD%"
set START=%TIME%

for /f "tokens=3" %%v in ('findstr /r /c:"^project(lima VERSION" CMakeLists.txt') do set "VERSION=%%v"
if "%VERSION%"=="" call :fail "Could not read the version from CMakeLists.txt" & exit /b 1

rem ---------------------------------------------------------------- Version
echo Current version in CMakeLists.txt: %VERSION%
if %DRY_RUN%==1 goto :version_chosen
if %UPLOAD_ONLY%==1 goto :version_chosen
set "NEW_VERSION="
set /p "NEW_VERSION=Version to release (Enter keeps %VERSION%): "
if "%NEW_VERSION%"=="" set "NEW_VERSION=%VERSION%"
echo %NEW_VERSION%| findstr /r /x "[0-9][0-9]*\.[0-9][0-9]*\.[0-9][0-9]*" > nul || (call :fail "'%NEW_VERSION%' is not a version like 1.4.0" & exit /b 1)
set TAGGED=
for /f %%t in ('git ls-remote --tags origin refs/tags/v%NEW_VERSION%') do set TAGGED=1
if defined TAGGED call :fail "v%NEW_VERSION% is already released" & exit /b 1
if "%NEW_VERSION%"=="%VERSION%" goto :version_chosen
set DIRTY=
for /f %%s in ('git status --porcelain') do set DIRTY=1
if defined DIRTY call :fail "Commit your changes before bumping the version" & exit /b 1
echo Bumping the version to %NEW_VERSION%, and committing and pushing that
powershell -NoProfile -Command "(Get-Content CMakeLists.txt) -replace '^project\(lima VERSION [0-9.]+','project(lima VERSION %NEW_VERSION%' | Set-Content -Encoding ascii CMakeLists.txt" || (call :fail "Could not update CMakeLists.txt" & exit /b 1)
git commit --quiet -m "Release v%NEW_VERSION%" CMakeLists.txt || (call :fail "Committing the version bump failed" & exit /b 1)
git push --quiet origin HEAD || (call :fail "Pushing the version bump failed" & exit /b 1)
set "VERSION=%NEW_VERSION%"
:version_chosen

for /f %%c in ('git rev-parse HEAD') do set "COMMIT=%%c"
set "OUT=%REPO%\distribution\out\%VERSION%"
if %UPLOAD_ONLY%==1 (
    if not exist "%OUT%\commit.txt" call :fail "Nothing to upload, %OUT% has no finished release. Run release.bat first" & exit /b 1
    set /p COMMIT=<"%OUT%\commit.txt"
)
echo ### LIMA %VERSION% from commit %COMMIT%

rem ---------------------------------------------------------------- Preflight
for /f %%s in ('git status --porcelain') do set DIRTY=1
if %DRY_RUN%==0 (
    if defined DIRTY if %UPLOAD_ONLY%==0 call :fail "There are uncommitted changes. Commit them, or use --dry-run" & exit /b 1
    git fetch --quiet origin || (call :fail "git fetch failed" & exit /b 1)
    set ON_REMOTE=
    for /f %%b in ('git branch -r --contains %COMMIT%') do set ON_REMOTE=1
    if not defined ON_REMOTE call :fail "Commit %COMMIT% is not on GitHub yet. Push it first" & exit /b 1
    set TAGGED=
    for /f %%t in ('git ls-remote --tags origin refs/tags/v%VERSION%') do set TAGGED=1
    if defined TAGGED call :fail "v%VERSION% is already released. Bump the version in CMakeLists.txt" & exit /b 1
) else (
    if defined DIRTY echo Dry run with uncommitted changes: Windows uses the working tree, Linux uses %COMMIT%
)

for /f "delims=" %%p in ('wsl -d %WSL_DISTRO% -- wslpath -a "%REPO%"') do set "REPO_WSL=%%p"
for /f "delims=" %%p in ('wsl -d %WSL_DISTRO% -- wslpath -a "%OUT%"') do set "OUT_WSL=%%p"
set "SCRIPTS_WSL=%REPO_WSL%/distribution/linux"
set CHECK_ARGS=
if %DRY_RUN%==0 set CHECK_ARGS=--upload
wsl -d %WSL_DISTRO% -- bash "%SCRIPTS_WSL%/check-wsl.sh" %CHECK_ARGS% || (call :fail "WSL is not set up" & exit /b 1)

if %UPLOAD_ONLY%==1 goto :release

if exist "%OUT%" rmdir /s /q "%OUT%"
mkdir "%OUT%"
> "%OUT%\commit.txt" echo %COMMIT%

rem ---------------------------------------------------------------- Windows
echo.
echo ### Building for Windows, CUDA architectures %CUDA_ARCHITECTURES%
call "%VS_DIR%\VC\Auxiliary\Build\vcvars64.bat" > nul
where cl > nul 2>&1 || (call :fail "Visual Studio environment not found in %VS_DIR%" & exit /b 1)
set "CMAKE=%VS_DIR%\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe"
set "NINJA=%VS_DIR%\Common7\IDE\CommonExtensions\Microsoft\CMake\Ninja\ninja.exe"
set "BUILD=%REPO%\build\release-windows"
"%CMAKE%" -S "%REPO%" -B "%BUILD%" -G Ninja -DCMAKE_MAKE_PROGRAM="%NINJA%" -DCMAKE_BUILD_TYPE=Release "-DCMAKE_CUDA_ARCHITECTURES=%CUDA_ARCHITECTURES%" || (call :fail "CMake configure failed" & exit /b 1)
"%CMAKE%" --build "%BUILD%" --target lima limaclitest || (call :fail "Windows build failed" & exit /b 1)

powershell -NoProfile -ExecutionPolicy Bypass -File "%REPO%\distribution\windows\package.ps1" -Executable "%BUILD%\code\LIMA\lima.exe" -OutDir "%OUT%" -Version %VERSION% -Commit %COMMIT% || (call :fail "Windows packaging failed" & exit /b 1)

echo.
echo ### Testing the packaged Windows build
"%BUILD%\code\LIMA_TESTS\limaclitest.exe" --lima "%OUT%\lima-%VERSION%-windows-x64\lima.exe" || (call :fail "The CLI tests failed on the packaged Windows build" & exit /b 1)

rem ---------------------------------------------------------------- Linux
echo.
echo ### Building for Linux in WSL
wsl -d %WSL_DISTRO% -- bash "%SCRIPTS_WSL%/build.sh" "%REPO_WSL%" %COMMIT% %VERSION% "%OUT_WSL%" || (call :fail "Linux build or smoke test failed" & exit /b 1)
echo.
echo ### Debian package
wsl -d %WSL_DISTRO% -- bash "%SCRIPTS_WSL%/package-deb.sh" %VERSION% "%OUT_WSL%" || (call :fail ".deb packaging or install test failed" & exit /b 1)
echo.
echo ### Arch package
wsl -d %WSL_DISTRO% -- bash "%SCRIPTS_WSL%/package-arch.sh" %VERSION% "%OUT_WSL%" || (call :fail "PKGBUILD generation or install test failed" & exit /b 1)

rem ---------------------------------------------------------------- Release
:release
for %%f in ("%OUT%\lima-%VERSION%-windows-x64.zip" "%OUT%\lima-%VERSION%-linux-x86_64.tar.gz" "%OUT%\lima_%VERSION%_amd64.deb" "%OUT%\PKGBUILD") do (
    if not exist %%f call :fail "Missing %%~nxf in %OUT%" & exit /b 1
)
powershell -NoProfile -Command "(Get-Content '%REPO%\distribution\release-notes.md') -replace '@VERSION@','%VERSION%' | Set-Content -Encoding utf8 '%OUT%\release-notes.md'"
set "ASSETS=%OUT_WSL%/lima-%VERSION%-windows-x64.zip %OUT_WSL%/lima-%VERSION%-linux-x86_64.tar.gz %OUT_WSL%/lima_%VERSION%_amd64.deb %OUT_WSL%/PKGBUILD"

if %DRY_RUN%==1 (
    echo.
    echo ### Dry run complete, nothing was uploaded. Release files are in %OUT%
    goto :done
)

echo.
echo ### Uploading draft release v%VERSION%
wsl -d %WSL_DISTRO% -- gh release create v%VERSION% --repo %GITHUB_REPO% --draft --target %COMMIT% --title "LIMA %VERSION%" --notes-file "%OUT_WSL%/release-notes.md" --generate-notes %ASSETS% || (call :fail "Uploading the release failed" & exit /b 1)
echo.
echo ### Draft release v%VERSION% is uploaded. Review and publish it at:
echo     https://github.com/%GITHUB_REPO%/releases

:done
echo Started %START%, finished %TIME%
popd
exit /b 0

:fail
echo.
echo RELEASE FAILED: %~1
popd
exit /b 1
