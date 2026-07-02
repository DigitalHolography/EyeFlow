<#
.SYNOPSIS
Build the EyeFlow Windows installer.

.DESCRIPTION
Reads the application version from pyproject.toml, builds a PyInstaller
one-dir bundle, generates an Inno Setup script, and compiles the final setup
executable into dist\installer.

.EXAMPLE
.\build_installer.ps1

.EXAMPLE
.\build_installer.ps1 -InnoSetupCompiler "C:\Program Files (x86)\Inno Setup 6\ISCC.exe"
#>

[CmdletBinding()]
param(
    [string]$InnoSetupCompiler = "",
    [string]$Python = "",
    [switch]$IncludePipelineExtras,
    [switch]$Console,
    [switch]$SkipClean
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"
$PSNativeCommandUseErrorActionPreference = $false

$RepoRoot = $PSScriptRoot
$PyprojectPath = Join-Path $RepoRoot "pyproject.toml"
$BuildRoot = Join-Path $RepoRoot "build\installer"
$PyInstallerWorkDir = Join-Path $BuildRoot "pyinstaller-work"
$PyInstallerSpecDir = Join-Path $BuildRoot "pyinstaller-spec"
$PyInstallerDistDir = Join-Path $BuildRoot "pyinstaller-dist"
$InstallerOutputDir = Join-Path $RepoRoot "dist\installer"
$GeneratedEntryPoint = Join-Path $BuildRoot "eyeflow_gui_entry.py"
$GeneratedInnoScript = Join-Path $BuildRoot "EyeFlow.iss"

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

    if (-not $targetFull.StartsWith($prefix, [System.StringComparison]::OrdinalIgnoreCase)) {
        throw "Refusing to operate outside the repository: $targetFull"
    }
}

function Read-ProjectMetadata {
    if (-not (Test-Path -LiteralPath $PyprojectPath -PathType Leaf)) {
        throw "pyproject.toml was not found at $PyprojectPath"
    }

    $content = Get-Content -LiteralPath $PyprojectPath -Raw -Encoding UTF8
    $name = "EyeFlow"
    $version = $null

    if ($content -match '(?m)^\s*name\s*=\s*"([^"]+)"\s*$') {
        $name = $Matches[1]
    }
    if ($content -match '(?m)^\s*version\s*=\s*"([^"]+)"\s*$') {
        $version = $Matches[1]
    }
    if (-not $version) {
        throw "Could not read [project].version from pyproject.toml"
    }

    return [pscustomobject]@{
        Name = $name
        Version = $version
    }
}

function ConvertTo-InnoQuotedValue {
    param([Parameter(Mandatory = $true)][string]$Value)
    return $Value.Replace('"', '""')
}

function ConvertTo-InnoVersionInfo {
    param([Parameter(Mandatory = $true)][string]$Version)

    if ($Version -notmatch '^\d+(\.\d+){0,3}$') {
        return $null
    }

    $parts = New-Object System.Collections.Generic.List[string]
    foreach ($part in $Version.Split(".")) {
        $parts.Add($part)
    }
    while ($parts.Count -lt 4) {
        $parts.Add("0")
    }
    return ($parts -join ".")
}

function ConvertTo-VersionDirName {
    param([Parameter(Mandatory = $true)][string]$Version)

    $safeVersion = ($Version -replace '[<>:"/\\|?*]+', '-').TrimEnd([char[]]" .")
    if ([string]::IsNullOrWhiteSpace($safeVersion)) {
        throw "Version cannot be converted to an installer directory name."
    }
    if ($safeVersion.StartsWith("v", [System.StringComparison]::OrdinalIgnoreCase)) {
        return $safeVersion
    }
    return "v$safeVersion"
}

function Join-AppNameAndVersion {
    param(
        [Parameter(Mandatory = $true)][string]$AppName,
        [Parameter(Mandatory = $true)][string]$AppVersion
    )

    if ([string]::IsNullOrWhiteSpace($AppVersion)) {
        return $AppName
    }
    return "$AppName $AppVersion"
}

function Resolve-InnoSetupCompiler {
    if ($InnoSetupCompiler) {
        if (-not (Test-Path -LiteralPath $InnoSetupCompiler -PathType Leaf)) {
            throw "Inno Setup compiler was not found: $InnoSetupCompiler"
        }
        return (Get-FullPath $InnoSetupCompiler)
    }

    $pathCommand = Get-Command "ISCC.exe" -ErrorAction SilentlyContinue
    if ($pathCommand) {
        return $pathCommand.Source
    }

    $candidates = @(
        "${env:ProgramFiles(x86)}\Inno Setup 6\ISCC.exe",
        "${env:ProgramFiles}\Inno Setup 6\ISCC.exe"
    )
    foreach ($candidate in $candidates) {
        if ($candidate -and (Test-Path -LiteralPath $candidate -PathType Leaf)) {
            return (Get-FullPath $candidate)
        }
    }

    throw "Inno Setup compiler (ISCC.exe) was not found. Install Inno Setup 6 or pass -InnoSetupCompiler."
}

function Resolve-PythonExe {
    if ($Python) {
        if (-not (Test-Path -LiteralPath $Python -PathType Leaf)) {
            throw "Python executable was not found: $Python"
        }
        return (Get-FullPath $Python)
    }

    $venvPython = Join-Path $RepoRoot ".venv\Scripts\python.exe"
    if (Test-Path -LiteralPath $venvPython -PathType Leaf) {
        return (Get-FullPath $venvPython)
    }

    $pathCommand = Get-Command "python.exe" -ErrorAction SilentlyContinue
    if ($pathCommand) {
        return $pathCommand.Source
    }

    throw "Python was not found. Install Python, run uv sync, or pass -Python."
}

function Invoke-PyInstaller {
    param([Parameter(Mandatory = $true)][string[]]$Arguments)

    $uv = Get-Command "uv.exe" -ErrorAction SilentlyContinue
    if (-not $uv) {
        $uv = Get-Command "uv" -ErrorAction SilentlyContinue
    }

    if ($uv -and -not $Python) {
        $uvArgs = @("run")
        if ($IncludePipelineExtras) {
            $uvArgs += @("--extra", "pipelines")
        }
        $uvArgs += @("--with", "pyinstaller", "python", "-m", "PyInstaller")
        $uvArgs += $Arguments
        & $uv.Source @uvArgs
        if ($LASTEXITCODE -ne 0) {
            throw "PyInstaller failed with exit code $LASTEXITCODE"
        }
        return
    }

    $pythonExe = Resolve-PythonExe
    & $pythonExe -m PyInstaller @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "PyInstaller failed with exit code $LASTEXITCODE"
    }
}

function Write-GeneratedEntryPoint {
    New-Item -ItemType Directory -Force -Path $BuildRoot | Out-Null
    @'
from launcher import main


if __name__ == "__main__":
    main()
'@ | Set-Content -LiteralPath $GeneratedEntryPoint -Encoding UTF8
}

function Write-InnoSetupScript {
    param(
        [Parameter(Mandatory = $true)][string]$AppName,
        [Parameter(Mandatory = $true)][string]$AppVersion,
        [Parameter(Mandatory = $true)][string]$AppDisplayName,
        [Parameter(Mandatory = $true)][string]$InstallRootName,
        [Parameter(Mandatory = $true)][string]$VersionDirName,
        [Parameter(Mandatory = $true)][string]$BundleDir
    )

    $appNameInno = ConvertTo-InnoQuotedValue $AppName
    $appVersionInno = ConvertTo-InnoQuotedValue $AppVersion
    $appDisplayNameInno = ConvertTo-InnoQuotedValue $AppDisplayName
    $installRootNameInno = ConvertTo-InnoQuotedValue $InstallRootName
    $versionDirNameInno = ConvertTo-InnoQuotedValue $VersionDirName
    $bundleDirInno = ConvertTo-InnoQuotedValue (Get-FullPath $BundleDir)
    $outputDirInno = ConvertTo-InnoQuotedValue (Get-FullPath $InstallerOutputDir)
    $licensePathInno = ConvertTo-InnoQuotedValue (Get-FullPath (Join-Path $RepoRoot "LICENSE"))
    $readmePathInno = ConvertTo-InnoQuotedValue (Get-FullPath (Join-Path $RepoRoot "README.md"))
    $thirdPartyNoticesPathInno = ConvertTo-InnoQuotedValue (Get-FullPath (Join-Path $RepoRoot "THIRD_PARTY_NOTICES"))
    $iconPathInno = ConvertTo-InnoQuotedValue (Get-FullPath (Join-Path $RepoRoot "EyeFlow.ico"))
    $pipelinesSourceDirInno = ConvertTo-InnoQuotedValue (Get-FullPath (Join-Path $RepoRoot "src\pipelines"))
    $setupBaseNameInno = ConvertTo-InnoQuotedValue ("$AppName-Setup-$AppVersion")
    $versionInfo = ConvertTo-InnoVersionInfo $AppVersion
    $versionInfoLine = ""
    if ($versionInfo) {
        $versionInfoLine = "VersionInfoVersion=$versionInfo"
    }

    @"
#define MyAppName "$appNameInno"
#define MyAppVersion "$appVersionInno"
#define MyAppDisplayName "$appDisplayNameInno"
#define MyInstallRootName "$installRootNameInno"
#define MyVersionDirName "$versionDirNameInno"
#define MyAppPublisher "EyeFlow"
#define MyAppExeName "$appDisplayNameInno.exe"
#define MyBundleDir "$bundleDirInno"
#define MyOutputDir "$outputDirInno"
#define MySetupBaseName "$setupBaseNameInno"
#define MyLicensePath "$licensePathInno"
#define MyReadmePath "$readmePathInno"
#define MyThirdPartyNoticesPath "$thirdPartyNoticesPathInno"
#define MyIconPath "$iconPathInno"
#define MyPipelinesSourceDir "$pipelinesSourceDirInno"

[Setup]
AppId={{72E10B9F-83C4-4F4C-89A5-F8BB598B0396}
AppName={#MyAppName}
AppVersion={#MyAppVersion}
AppPublisher={#MyAppPublisher}
AppVerName={#MyAppDisplayName}
DefaultDirName={localappdata}\Programs\{#MyInstallRootName}\{#MyVersionDirName}
DefaultGroupName={#MyAppName}
DisableProgramGroupPage=yes
UsePreviousAppDir=no
LicenseFile={#MyLicensePath}
OutputDir={#MyOutputDir}
OutputBaseFilename={#MySetupBaseName}
SetupIconFile={#MyIconPath}
UninstallDisplayIcon={app}\{#MyAppExeName}
Compression=lzma2
SolidCompression=yes
WizardStyle=modern
PrivilegesRequired=lowest
ArchitecturesAllowed=x64compatible
$versionInfoLine

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Tasks]
Name: "desktopicon"; Description: "Create a desktop shortcut"; GroupDescription: "Additional icons:"; Flags: unchecked

[Files]
Source: "{#MyBundleDir}\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs
Source: "{#MyLicensePath}"; DestDir: "{app}"; DestName: "LICENSE"; Flags: ignoreversion
Source: "{#MyReadmePath}"; DestDir: "{app}"; DestName: "README.md"; Flags: ignoreversion
Source: "{#MyThirdPartyNoticesPath}"; DestDir: "{app}"; DestName: "THIRD_PARTY_NOTICES"; Flags: ignoreversion
Source: "{#MyPipelinesSourceDir}\*"; DestDir: "{app}\pipelines"; Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
Name: "{group}\{#MyAppDisplayName}"; Filename: "{app}\{#MyAppExeName}"
Name: "{group}\Uninstall {#MyAppDisplayName}"; Filename: "{uninstallexe}"
Name: "{autodesktop}\{#MyAppDisplayName}"; Filename: "{app}\{#MyAppExeName}"; Tasks: desktopicon

[Run]
Filename: "{app}\{#MyAppExeName}"; Description: "Launch {#MyAppDisplayName}"; Flags: nowait postinstall skipifsilent
"@ | Set-Content -LiteralPath $GeneratedInnoScript -Encoding UTF8
}

$metadata = Read-ProjectMetadata
$appName = $metadata.Name
$appVersion = $metadata.Version
$appDisplayName = Join-AppNameAndVersion -AppName $appName -AppVersion $appVersion
$installRootName = "eyeflow-python"
$versionDirName = ConvertTo-VersionDirName -Version $appVersion
$bundleDir = Join-Path $PyInstallerDistDir $appDisplayName
$bundleExe = Join-Path $bundleDir "$appDisplayName.exe"

Write-Host "Building $appDisplayName installer..."

if (-not $SkipClean) {
    Assert-ChildPath -BasePath $RepoRoot -TargetPath $BuildRoot
    Assert-ChildPath -BasePath $RepoRoot -TargetPath $InstallerOutputDir
    Remove-Item -LiteralPath $BuildRoot -Recurse -Force -ErrorAction SilentlyContinue
    Remove-Item -LiteralPath $InstallerOutputDir -Recurse -Force -ErrorAction SilentlyContinue
}

New-Item -ItemType Directory -Force -Path $BuildRoot | Out-Null
New-Item -ItemType Directory -Force -Path $InstallerOutputDir | Out-Null
Write-GeneratedEntryPoint

$pyInstallerArgs = @(
    "--noconfirm",
    "--clean",
    "--onedir",
    "--name", $appDisplayName,
    "--distpath", $PyInstallerDistDir,
    "--workpath", $PyInstallerWorkDir,
    "--specpath", $PyInstallerSpecDir,
    "--paths", (Join-Path $RepoRoot "src"),
    "--icon", (Join-Path $RepoRoot "EyeFlow.ico"),
    "--add-data", "$((Join-Path $RepoRoot "EyeFlow_logo.png"));.",
    "--add-data", "$((Join-Path $RepoRoot "default_settings.json"));.",
    "--add-data", "$($PyprojectPath);.",
    "--hidden-import", "eye_flow",
    "--collect-submodules", "pipelines",
    "--collect-all", "tkinterdnd2",
    "--collect-all", "sv_ttk"
)

if ($Console) {
    $pyInstallerArgs += "--console"
} else {
    $pyInstallerArgs += "--windowed"
}

$pyInstallerArgs += $GeneratedEntryPoint

Invoke-PyInstaller -Arguments $pyInstallerArgs

if (-not (Test-Path -LiteralPath $bundleExe -PathType Leaf)) {
    throw "PyInstaller did not produce the expected executable: $bundleExe"
}

$iscc = Resolve-InnoSetupCompiler
Write-InnoSetupScript `
    -AppName $appName `
    -AppVersion $appVersion `
    -AppDisplayName $appDisplayName `
    -InstallRootName $installRootName `
    -VersionDirName $versionDirName `
    -BundleDir $bundleDir

& $iscc $GeneratedInnoScript
if ($LASTEXITCODE -ne 0) {
    throw "Inno Setup failed with exit code $LASTEXITCODE"
}

$installerPath = Join-Path $InstallerOutputDir "$appName-Setup-$appVersion.exe"
if (-not (Test-Path -LiteralPath $installerPath -PathType Leaf)) {
    throw "Inno Setup did not produce the expected installer: $installerPath"
}

Write-Host "Installer created: $installerPath"
