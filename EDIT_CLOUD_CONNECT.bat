@echo off
setlocal
cd /d "%~dp0"

set "RDIR="
for /d %%D in ("C:\Program Files\R\R-*") do set "RDIR=%%D"

if not defined RDIR (
  echo R installation was not found under C:\Program Files\R\
  pause
  exit /b 1
)

echo Using R from:
echo %RDIR%
echo.

"%RDIR%\bin\Rscript.exe" "scripts\prepare_connect_cloud.R"
if errorlevel 1 (
  echo.
  echo Failed to create manifest.json.
  pause
  exit /b 1
)

start "" "manifest.json"
exit /b 0
