@echo off
setlocal
cd /d "%~dp0"
python omega_v2.py
if errorlevel 1 pause