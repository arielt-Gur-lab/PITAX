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

echo [1/5] Taxonomy logic smoke tests
"%RDIR%\bin\Rscript.exe" "tests\taxonomy_logic_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [2/5] QC ambiguous-peak smoke tests
"%RDIR%\bin\Rscript.exe" "tests\qc_peak_flags_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [3/5] Manual curation / undo-redo smoke tests
"%RDIR%\bin\Rscript.exe" "tests\manual_curation_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [4/5] Stage 1 AB1 evidence helper tests
"%RDIR%\bin\Rscript.exe" "tests\ab1_evidence_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [5/5] Workflow order and QC naming smoke tests
"%RDIR%\bin\Rscript.exe" "tests\workflow_order_smoke.R"
if errorlevel 1 goto :failed

echo.
echo All PITAX v3.0.0-alpha.4 tests passed.
pause
exit /b 0

:failed
echo.
echo One or more tests FAILED.
pause
exit /b 1
