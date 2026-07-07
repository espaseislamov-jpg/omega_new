@echo off
setlocal
cd /d "%~dp0\.."
py -3 -m pip install --upgrade pip
py -3 -m pip install -r requirements.txt
py -3 scripts\build_omega_v2.py
endlocal
