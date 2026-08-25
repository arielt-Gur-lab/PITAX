@echo off
setlocal
cd /d "%~dp0"

set "PITAX_VERSION="
set /p PITAX_VERSION=<VERSION.txt
if not defined PITAX_VERSION set "PITAX_VERSION=unknown"

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

echo [preflight] Parse all R sources
"%RDIR%\bin\Rscript.exe" "tests\preflight\parse_all_r_sources.R"
if errorlevel 1 goto :failed

echo.
echo [1/16] Taxonomy logic smoke tests
"%RDIR%\bin\Rscript.exe" "tests\unit\taxonomy_logic_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [2/16] QC ambiguous-peak smoke tests
"%RDIR%\bin\Rscript.exe" "tests\unit\qc_peak_flags_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [3/16] Manual curation / undo-redo smoke tests
"%RDIR%\bin\Rscript.exe" "tests\unit\manual_curation_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [4/16] Stage 1 AB1 evidence helper tests
"%RDIR%\bin\Rscript.exe" "tests\unit\ab1_evidence_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [5/16] Workflow order and QC naming smoke tests
"%RDIR%\bin\Rscript.exe" "tests\contracts\workflow_order_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [6/16] Stage 2 architecture and migration tests
"%RDIR%\bin\Rscript.exe" "tests\unit\stage2_architecture_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [7/16] Stage 2 closed-gate Shiny integration tests
"%RDIR%\bin\Rscript.exe" "tests\contracts\stage2_app_contract_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [8/16] Stage 3 consensus algorithm tests
"%RDIR%\bin\Rscript.exe" "tests\unit\stage3_consensus_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [9/16] Stage 3 controlled AB1 pair tests
"%RDIR%\bin\Rscript.exe" "tests\integration\stage3_ab1_fixture_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [10/16] Stage 3 Shiny integration and startup tests
"%RDIR%\bin\Rscript.exe" "tests\integration\stage3_app_contract_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [11/16] Stage 4 multi-locus profile tests
"%RDIR%\bin\Rscript.exe" "tests\unit\stage4_multilocus_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [12/16] Stage 4 Shiny integration contract
"%RDIR%\bin\Rscript.exe" "tests\contracts\stage4_app_contract_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [13/16] Global DataTables alignment regression contract
"%RDIR%\bin\Rscript.exe" "tests\contracts\datatables_alignment_contract_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [14/16] Organized project structure contract
"%RDIR%\bin\Rscript.exe" "tests\contracts\project_structure_contract_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [15/16] Alpha 10 assay and schema 6 foundation tests
"%RDIR%\bin\Rscript.exe" "tests\unit\alpha10_assay_schema_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [16/16] Reflection follow-up correctness tests
"%RDIR%\bin\Rscript.exe" "tests\unit\reflection_followup_smoke.R"
if errorlevel 1 goto :failed

echo.
echo All PITAX v%PITAX_VERSION% tests passed.
pause
exit /b 0

:failed
echo.
echo One or more tests FAILED.
pause
exit /b 1
