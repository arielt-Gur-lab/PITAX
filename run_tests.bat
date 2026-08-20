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

echo [1/10] Taxonomy logic smoke tests
"%RDIR%\bin\Rscript.exe" "tests\taxonomy_logic_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [2/10] QC ambiguous-peak smoke tests
"%RDIR%\bin\Rscript.exe" "tests\qc_peak_flags_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [3/10] Manual curation / undo-redo smoke tests
"%RDIR%\bin\Rscript.exe" "tests\manual_curation_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [4/10] Stage 1 AB1 evidence helper tests
"%RDIR%\bin\Rscript.exe" "tests\ab1_evidence_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [5/10] Workflow order and QC naming smoke tests
"%RDIR%\bin\Rscript.exe" "tests\workflow_order_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [6/10] Stage 2 architecture and migration tests
"%RDIR%\bin\Rscript.exe" "tests\stage2_architecture_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [7/10] Stage 2 closed-gate Shiny integration tests
"%RDIR%\bin\Rscript.exe" "tests\stage2_app_contract_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [8/10] Stage 3 consensus algorithm tests
"%RDIR%\bin\Rscript.exe" "tests\stage3_consensus_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [9/10] Stage 3 controlled AB1 pair tests
"%RDIR%\bin\Rscript.exe" "tests\stage3_ab1_fixture_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [10/10] Stage 3 Shiny integration and startup tests
"%RDIR%\bin\Rscript.exe" "tests\stage3_app_contract_smoke.R"
if errorlevel 1 goto :failed

echo.
echo All PITAX v3.0.0-alpha.8.2 tests passed.
pause
exit /b 0

:failed
echo.
echo One or more tests FAILED.
pause
exit /b 1
