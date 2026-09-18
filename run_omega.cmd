@echo off
setlocal
cd /d "%~dp0"
if exist ".venv\Scripts\python.exe" (
    ".venv\Scripts\python.exe" omega_v2.py
) else (
    python omega_v2.py
)
if errorlevel 1 pause
