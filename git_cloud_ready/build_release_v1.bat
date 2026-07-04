@echo off
setlocal

set APP_NAME=Omega_v1
set SCRIPT=New_idea.py
set RELEASE_DIR=release\%APP_NAME%
set ZIP_NAME=release\%APP_NAME%_release.zip

if exist build rmdir /s /q build
if exist dist\%APP_NAME% rmdir /s /q dist\%APP_NAME%
if exist %RELEASE_DIR% rmdir /s /q %RELEASE_DIR%
if exist %ZIP_NAME% del /f /q %ZIP_NAME%

py -3 -m PyInstaller ^
  --noconfirm ^
  --clean ^
  --onedir ^
  --windowed ^
  --name %APP_NAME% ^
  --add-data "reference_targets_reverted_c22fixed.json;." ^
  --add-data "chromatogram_gui_settings.json;." ^
  --collect-all pyopenms ^
  --collect-all lmfit ^
  --hidden-import=tkinter ^
  --hidden-import=matplotlib.backends.backend_tkagg ^
  %SCRIPT%

mkdir %RELEASE_DIR%
xcopy /e /i /y dist\%APP_NAME% %RELEASE_DIR% >nul
copy /y README_v1.txt %RELEASE_DIR%\README_v1.txt >nul
copy /y reference_targets_reverted_c22fixed.json %RELEASE_DIR%\reference_targets_reverted_c22fixed.json >nul
copy /y chromatogram_gui_settings.json %RELEASE_DIR%\chromatogram_gui_settings.json >nul

py -3 -c "from pathlib import Path; import zipfile; root = Path(r'%RELEASE_DIR%'); zip_path = Path(r'%ZIP_NAME%'); zip_path.parent.mkdir(parents=True, exist_ok=True); \
with zipfile.ZipFile(zip_path, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as zf: \
    [zf.write(path, path.relative_to(root.parent)) for path in root.rglob('*') if path.is_file()]"

endlocal
