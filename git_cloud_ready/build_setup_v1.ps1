$ErrorActionPreference = "Stop"

$root = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $root

$appName = "Omega_v1"
$distDir = Join-Path $root "dist\$appName"
$buildDir = Join-Path $root "build"
$releaseDir = Join-Path $root "release\$appName"
$payloadZip = Join-Path $root "release\${appName}_payload.zip"
$setupExe = Join-Path $root "release\${appName}_setup.exe"
$setupSource = Join-Path $root "SetupInstaller.cs"
$cscPath = Join-Path $env:WINDIR "Microsoft.NET\Framework64\v4.0.30319\csc.exe"
$legacyInstallPs1 = Join-Path $root "release\install_release.ps1"
$legacyInstallVbs = Join-Path $root "release\install_release.vbs"
$legacySed = Join-Path $root "release\${appName}_setup.sed"
$legacyDdf = Join-Path $root "release\~${appName}_setup.DDF"

if (Test-Path $buildDir) {
    Remove-Item $buildDir -Recurse -Force
}
if (Test-Path $distDir) {
    Remove-Item $distDir -Recurse -Force
}
if (Test-Path $releaseDir) {
    Remove-Item $releaseDir -Recurse -Force
}
if (Test-Path $payloadZip) {
    Remove-Item $payloadZip -Force
}
if (Test-Path $setupExe) {
    Remove-Item $setupExe -Force
}
if (Test-Path $legacyInstallPs1) {
    Remove-Item $legacyInstallPs1 -Force
}
if (Test-Path $legacyInstallVbs) {
    Remove-Item $legacyInstallVbs -Force
}
if (Test-Path $legacySed) {
    Remove-Item $legacySed -Force
}
if (Test-Path $legacyDdf) {
    Remove-Item $legacyDdf -Force
}
py -3 -m PyInstaller --noconfirm --clean Omega_v1.spec
if ($LASTEXITCODE -ne 0) {
    throw "PyInstaller build failed."
}

New-Item -ItemType Directory -Path $releaseDir -Force | Out-Null
Copy-Item (Join-Path $distDir "*") $releaseDir -Recurse -Force
Copy-Item (Join-Path $root "README_v1.txt") (Join-Path $releaseDir "README_v1.txt") -Force
Copy-Item (Join-Path $root "reference_targets_reverted_c22fixed.json") (Join-Path $releaseDir "reference_targets_reverted_c22fixed.json") -Force
Copy-Item (Join-Path $root "chromatogram_gui_settings.json") (Join-Path $releaseDir "chromatogram_gui_settings.json") -Force
Copy-Item (Join-Path $root "uninstall_release.ps1") (Join-Path $releaseDir "uninstall_release.ps1") -Force

Compress-Archive -Path (Join-Path $releaseDir "*") -DestinationPath $payloadZip -CompressionLevel Optimal -Force
& $cscPath /nologo /target:winexe /platform:x64 /optimize+ /out:$setupExe /resource:$payloadZip,Omega_v1_payload.zip /r:System.dll /r:System.Windows.Forms.dll /r:System.Drawing.dll /r:System.IO.Compression.dll /r:System.IO.Compression.FileSystem.dll $setupSource
if ($LASTEXITCODE -ne 0) {
    throw "Setup build failed."
}

Write-Host "Setup ready:" $setupExe
