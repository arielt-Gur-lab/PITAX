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

echo [1/3] Taxonomy logic smoke tests
"%RDIR%\bin\Rscript.exe" "tests\taxonomy_logic_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [2/3] QC ambiguous-peak smoke tests
"%RDIR%\bin\Rscript.exe" "tests\qc_peak_flags_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [3/3] Manual curation / undo-redo smoke tests
"%RDIR%\bin\Rscript.exe" "tests\manual_curation_smoke.R"
if errorlevel 1 goto :failed

echo.
echo All v2.14.0 smoke tests passed.
pause
exit /b 0

:failed
echo.
echo One or more smoke tests FAILED.
pause
exit /b 1
