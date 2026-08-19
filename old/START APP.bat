@echo off
cd /d "%~dp0"

for /d %%D in ("C:\Program Files\R\R-*") do (
    set "RDIR=%%D"
)

if not defined RDIR (
    echo R installation was not found.
    pause
    exit /b 1
)

echo Using R from:
echo %RDIR%

"%RDIR%\bin\Rscript.exe" -e "shiny::runApp('.', launch.browser=TRUE)"

pause