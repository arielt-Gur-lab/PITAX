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

echo [1/15] Taxonomy logic smoke tests
"%RDIR%\bin\Rscript.exe" "tests\unit\taxonomy_logic_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [2/15] QC ambiguous-peak smoke tests
"%RDIR%\bin\Rscript.exe" "tests\unit\qc_peak_flags_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [3/15] Manual curation / undo-redo smoke tests
"%RDIR%\bin\Rscript.exe" "tests\unit\manual_curation_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [4/15] Stage 1 AB1 evidence helper tests
"%RDIR%\bin\Rscript.exe" "tests\unit\ab1_evidence_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [5/15] Workflow order and QC naming smoke tests
"%RDIR%\bin\Rscript.exe" "tests\contracts\workflow_order_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [6/15] Stage 2 architecture and migration tests
"%RDIR%\bin\Rscript.exe" "tests\unit\stage2_architecture_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [7/15] Stage 2 closed-gate Shiny integration tests
"%RDIR%\bin\Rscript.exe" "tests\contracts\stage2_app_contract_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [8/15] Stage 3 consensus algorithm tests
"%RDIR%\bin\Rscript.exe" "tests\unit\stage3_consensus_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [9/15] Stage 3 controlled AB1 pair tests
"%RDIR%\bin\Rscript.exe" "tests\integration\stage3_ab1_fixture_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [10/15] Stage 3 Shiny integration and startup tests
"%RDIR%\bin\Rscript.exe" "tests\integration\stage3_app_contract_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [11/15] Stage 4 multi-locus profile tests
"%RDIR%\bin\Rscript.exe" "tests\unit\stage4_multilocus_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [12/15] Stage 4 Shiny integration contract
"%RDIR%\bin\Rscript.exe" "tests\contracts\stage4_app_contract_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [13/15] Global DataTables alignment regression contract
"%RDIR%\bin\Rscript.exe" "tests\contracts\datatables_alignment_contract_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [14/15] Organized project structure contract
"%RDIR%\bin\Rscript.exe" "tests\contracts\project_structure_contract_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [15/15] Alpha 10 assay and schema 6 foundation tests
"%RDIR%\bin\Rscript.exe" "tests\unit\alpha10_assay_schema_smoke.R"
if errorlevel 1 goto :failed

echo.
echo All PITAX v3.0.0-alpha.10.3 tests passed.
pause
exit /b 0

:failed
echo.
echo One or more tests FAILED.
pause
exit /b 1
