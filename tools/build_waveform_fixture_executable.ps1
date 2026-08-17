<#
.SYNOPSIS
Build the standalone EyeFlow Waveform Fixture Extractor for Windows.

.DESCRIPTION
Uses PyInstaller to create a single, windowed executable and writes a matching
SHA-256 checksum file. The selected Python environment must already contain
numpy, h5py, and pyinstaller.
#>

[CmdletBinding()]
param(
    [string]$Python = "",
    [switch]$SkipClean
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"
$PSNativeCommandUseErrorActionPreference = $false

$RepoRoot = Split-Path -Parent $PSScriptRoot
$BuildRoot = Join-Path $RepoRoot "build\waveform-fixture-extractor"
$DistRoot = Join-Path $RepoRoot "dist\waveform-fixture-extractor"
$ExecutableName = "EyeFlowWaveformFixtureExtractor.exe"
$ExecutablePath = Join-Path $DistRoot $ExecutableName
$ChecksumPath = "$ExecutablePath.sha256"

function Get-FullPath {
    param([Parameter(Mandatory = $true)][string]$Path)
    return [System.IO.Path]::GetFullPath($Path)
}

function Assert-ChildPath {
    param(
        [Parameter(Mandatory = $true)][string]$BasePath,
        [Parameter(Mandatory = $true)][string]$TargetPath
    )

    $baseFull = (Get-FullPath $BasePath).TrimEnd(
        [System.IO.Path]::DirectorySeparatorChar,
        [System.IO.Path]::AltDirectorySeparatorChar
    )
    $targetFull = Get-FullPath $TargetPath
    $prefix = $baseFull + [System.IO.Path]::DirectorySeparatorChar
    if (-not $targetFull.StartsWith(
        $prefix,
        [System.StringComparison]::OrdinalIgnoreCase
    )) {
        throw "Refusing to operate outside the repository: $targetFull"
    }
}

function Resolve-Python {
    if ($Python) {
        if (-not (Test-Path -LiteralPath $Python -PathType Leaf)) {
            throw "Python executable was not found: $Python"
        }
        return (Get-FullPath $Python)
    }

    $command = Get-Command "python.exe" -ErrorAction SilentlyContinue
    if (-not $command) {
        throw "Python was not found. Pass -Python with an explicit path."
    }
    return $command.Source
}

Assert-ChildPath -BasePath $RepoRoot -TargetPath $BuildRoot
Assert-ChildPath -BasePath $RepoRoot -TargetPath $DistRoot
if (-not $SkipClean) {
    Remove-Item -LiteralPath $BuildRoot -Recurse -Force -ErrorAction SilentlyContinue
    Remove-Item -LiteralPath $DistRoot -Recurse -Force -ErrorAction SilentlyContinue
}
New-Item -ItemType Directory -Force -Path $BuildRoot | Out-Null
New-Item -ItemType Directory -Force -Path $DistRoot | Out-Null

$pythonExe = Resolve-Python
$entryPoint = Join-Path $PSScriptRoot "extract_waveform_fixture_gui.py"
$iconPath = Join-Path $RepoRoot "EyeFlow.ico"
$pyInstallerArgs = @(
    "-m", "PyInstaller",
    "--noconfirm",
    "--clean",
    "--onefile",
    "--windowed",
    "--noupx",
    "--name", "EyeFlowWaveformFixtureExtractor",
    "--distpath", $DistRoot,
    "--workpath", (Join-Path $BuildRoot "work"),
    "--specpath", (Join-Path $BuildRoot "spec"),
    "--paths", $PSScriptRoot,
    "--icon", $iconPath,
    "--hidden-import", "tkinter.filedialog",
    "--hidden-import", "tkinter.messagebox",
    $entryPoint
)

& $pythonExe @pyInstallerArgs
if ($LASTEXITCODE -ne 0) {
    throw "PyInstaller failed with exit code $LASTEXITCODE."
}
if (-not (Test-Path -LiteralPath $ExecutablePath -PathType Leaf)) {
    throw "Expected executable was not created: $ExecutablePath"
}

$hash = (Get-FileHash -LiteralPath $ExecutablePath -Algorithm SHA256).Hash.ToLowerInvariant()
"$hash *$ExecutableName" | Set-Content -LiteralPath $ChecksumPath -Encoding ASCII

Write-Host "Executable: $ExecutablePath"
Write-Host "SHA-256:   $hash"
Write-Host "Checksum:  $ChecksumPath"
