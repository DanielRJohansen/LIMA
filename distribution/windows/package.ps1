# Stages and zips the Windows release. Called by distribution/release.bat from a Visual Studio developer
# environment (needs dumpbin and VCToolsRedistDir).
#
# Layout of the zip:
#   lima-<version>-windows-x64/
#     lima.exe, cufft DLL, MSVC runtime DLLs
#     resources/
#     LICENSE.txt, README.txt
param(
    [Parameter(Mandatory)][string]$Executable,
    [Parameter(Mandatory)][string]$OutDir,
    [Parameter(Mandatory)][string]$Version,
    [Parameter(Mandatory)][string]$Commit
)
$ErrorActionPreference = 'Stop'
# Windows' own bsdtar, which also writes zips. GNU tar from Git for Windows may come first on PATH
$tar = Join-Path $env:SystemRoot 'System32\tar.exe'

$repo = (Resolve-Path "$PSScriptRoot\..\..").Path
$name = "lima-$Version-windows-x64"
$stage = Join-Path $OutDir $name
if (Test-Path $stage) { Remove-Item -Recurse -Force $stage }
New-Item -ItemType Directory -Force $stage | Out-Null

Copy-Item $Executable $stage

# Only git-tracked resources, from the commit being released (no caches or local experiments)
$resourcesTar = Join-Path $OutDir 'resources.tar'
git -C $repo archive --format=tar -o $resourcesTar $Commit resources
if ($LASTEXITCODE -ne 0) { throw 'git archive failed' }
& $tar -xf $resourcesTar -C $stage
if ($LASTEXITCODE -ne 0) { throw 'Extracting resources failed' }
Remove-Item $resourcesTar

function Get-Dependents([string]$binary) {
    $lines = & dumpbin /nologo /dependents $binary
    if ($LASTEXITCODE -ne 0) { throw "dumpbin failed on $binary" }
    return $lines | ForEach-Object { $_.Trim() } | Where-Object { $_ -match '^[\w.-]+\.dll$' }	# Bare names only, not the 'Dump of file' header
}

# cuFFT has no static library on Windows, so its DLL ships next to lima.exe (redistributable per the CUDA EULA)
$cudaPath = $env:CUDA_PATH
if (-not $cudaPath) { throw 'CUDA_PATH is not set' }
foreach ($dll in Get-Dependents (Join-Path $stage 'lima.exe') | Where-Object { $_ -like 'cufft*' }) {
    $source = Get-ChildItem -Path "$cudaPath\bin" -Recurse -Filter $dll | Select-Object -First 1
    if (-not $source) { throw "$dll not found under $cudaPath\bin" }
    Copy-Item $source.FullName $stage
}

# App-local MSVC runtime, so users dont need to install the Visual C++ redistributable
$redist = Get-ChildItem -Path "$env:VCToolsRedistDir\x64" -Directory -Filter 'Microsoft.VC*.CRT' | Select-Object -First 1
if (-not $redist) { throw "MSVC runtime not found in $env:VCToolsRedistDir\x64" }
Copy-Item "$($redist.FullName)\*.dll" $stage

Copy-Item (Join-Path $repo 'LICENSE.txt') $stage
Copy-Item (Join-Path $repo 'THIRD_PARTY_NOTICES.txt') $stage

# The citations for the bundled force fields and lipids must ship with every copy
$readme = git -C $repo show "${Commit}:README.md"
$citationsStart = ($readme | Select-String -Pattern '^## LIMA would not be possible' | Select-Object -First 1).LineNumber
if (-not $citationsStart) { throw 'The citations section was not found in README.md' }
$readme[($citationsStart - 1)..($readme.Count - 1)] | Set-Content (Join-Path $stage 'CITATIONS.md')
@"
LIMA $Version for Windows (x64)

Requirements
  An NVIDIA GPU of the RTX 40-series or newer (or H100, B200 and similar),
  with a recent NVIDIA driver (CUDA 13 capable).

Getting started
  Unzip anywhere and keep the files together; lima.exe finds the resources
  folder next to it. To run lima from any terminal, add this folder to your PATH.
  Run 'lima --help' for the available commands.

  Optional: run 'lima setregistry' once to open .gro files with LIMA by
  double-clicking. This only affects the current user and needs no admin rights.

License
  Free for small companies, noncommercial use and academia, and free to evaluate
  for 31 days, see LICENSE.txt. Other use requires a commercial license:
  daniel@lima-dynamics.com

https://github.com/DanielRJohansen/LIMA
"@ | Set-Content -Encoding ascii (Join-Path $stage 'README.txt')

# Every DLL that lima.exe or a shipped DLL loads must either ship with it or be part of Windows
$windowsDlls = @('KERNEL32.dll', 'USER32.dll', 'GDI32.dll', 'SHELL32.dll', 'ADVAPI32.dll', 'OPENGL32.dll', 'IMM32.dll',
    'ole32.dll', 'OLEAUT32.dll', 'WS2_32.dll', 'VERSION.dll', 'SETUPAPI.dll', 'dbghelp.dll', 'bcrypt.dll', 'CRYPT32.dll',
    'RPCRT4.dll', 'WINMM.dll', 'COMDLG32.dll', 'SHLWAPI.dll', 'dwmapi.dll', 'UxTheme.dll', 'ntdll.dll')
$shipped = Get-ChildItem $stage -Filter *.dll | ForEach-Object { $_.Name }
$missing = @()
foreach ($binary in @(Join-Path $stage 'lima.exe') + (Get-ChildItem $stage -Filter *.dll | ForEach-Object { $_.FullName })) {
    foreach ($dll in Get-Dependents $binary) {
        if ($dll -like 'api-ms-win-*' -or $dll -like 'ext-ms-win-*') { continue }
        if ($windowsDlls -contains $dll -or $shipped -contains $dll) { continue }
        $missing += "$dll (needed by $(Split-Path -Leaf $binary))"
    }
}
if ($missing) { throw "The release would not run on a clean machine, these DLLs are not shipped:`n  $($missing -join "`n  ")" }

$zip = Join-Path $OutDir "$name.zip"
if (Test-Path $zip) { Remove-Item $zip }
& $tar -a -c -f $zip -C $OutDir $name
if ($LASTEXITCODE -ne 0) { throw 'Creating the zip failed' }
Write-Host ("Packaged {0} ({1:N0} MB)" -f $zip, ((Get-Item $zip).Length / 1MB))
