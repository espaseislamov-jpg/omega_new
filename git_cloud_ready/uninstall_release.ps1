$ErrorActionPreference = "SilentlyContinue"

$appName = "Omega v1"
$installDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$desktopShortcut = Join-Path ([Environment]::GetFolderPath("Desktop")) "$appName.lnk"
$startMenuDir = Join-Path $env:APPDATA "Microsoft\Windows\Start Menu\Programs\$appName"

if (Test-Path $desktopShortcut) {
    Remove-Item $desktopShortcut -Force
}
if (Test-Path $startMenuDir) {
    Remove-Item $startMenuDir -Recurse -Force
}

$parentDir = Split-Path -Parent $installDir
if (Test-Path $installDir) {
    Remove-Item $installDir -Recurse -Force
}

if ((Test-Path $parentDir) -and -not (Get-ChildItem $parentDir -Force | Select-Object -First 1)) {
    Remove-Item $parentDir -Force
}
