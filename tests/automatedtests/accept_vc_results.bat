@echo off
setlocal

py -3 --version >nul 2>nul
if %errorlevel% equ 0 (
    py -3 "%~dp0accept_vc_results.py" %*
) else (
    python "%~dp0accept_vc_results.py" %*
)

set "accept_result=%errorlevel%"
if not "%accept_result%"=="0" (
    echo.
    echo Failed to accept VC results.
) else (
    echo.
    echo VC targets updated successfully.
)

if not "%~1"=="" exit /b %accept_result%
pause
exit /b %accept_result%
