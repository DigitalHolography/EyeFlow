<#
.SYNOPSIS
Prepare locked build dependencies and create the EyeFlow Windows installer.
.EXAMPLE
powershell ./build_installer.ps1
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

$builder = Join-Path $PSScriptRoot "installer\build_installer.py"
$builderArgs = @()
if ($InnoSetupCompiler) { $builderArgs += @("--inno-setup-compiler", $InnoSetupCompiler) }
if ($IncludePipelineExtras) { $builderArgs += "--include-pipeline-extras" }
if ($Console) { $builderArgs += "--console" }
if ($SkipClean) { $builderArgs += "--skip-clean" }

if (-not $Python) {
    $uv = Get-Command "uv" -ErrorAction SilentlyContinue
    if (-not $uv) {
        throw "Install uv, or pass -Python with the gpu and installer dependencies already installed."
    }
    # Keep build dependencies outside the cleaned staging directory.
    $environment = Join-Path $PSScriptRoot "build\installer-env"
    $previousEnvironment = $env:UV_PROJECT_ENVIRONMENT
    $env:UV_PROJECT_ENVIRONMENT = $environment
    try {
        $syncArgs = @("sync", "--project", $PSScriptRoot, "--locked", "--no-default-groups", "--extra", "gpu", "--extra", "installer")
        if ($IncludePipelineExtras) { $syncArgs += @("--extra", "pipelines") }
        & $uv.Source @syncArgs
        if ($LASTEXITCODE -ne 0) { throw "Build dependency preparation failed with exit code $LASTEXITCODE" }
    } finally {
        $env:UV_PROJECT_ENVIRONMENT = $previousEnvironment
    }
    $Python = Join-Path $environment "Scripts\python.exe"
}

& $Python $builder @builderArgs
if ($LASTEXITCODE -ne 0) { throw "Installer build failed with exit code $LASTEXITCODE" }
