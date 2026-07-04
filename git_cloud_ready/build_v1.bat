@echo off
setlocal

set APP_NAME=Omega_v1
set SCRIPT=New_idea.py

py -3 -m PyInstaller ^
  --noconfirm ^
  --clean ^
  --onefile ^
  --windowed ^
  --name %APP_NAME% ^
  --add-data "reference_targets_reverted_c22fixed.json;." ^
  --add-data "chromatogram_gui_settings.json;." ^
  --collect-all pyopenms ^
  --collect-all lmfit ^
  --collect-submodules matplotlib.backends ^
  --hidden-import=tkinter ^
  %SCRIPT%

endlocal
