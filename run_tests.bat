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
echo [1/19] Taxonomy logic smoke tests
"%RDIR%\bin\Rscript.exe" "tests\unit\taxonomy_logic_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [2/19] QC ambiguous-peak smoke tests
"%RDIR%\bin\Rscript.exe" "tests\unit\qc_peak_flags_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [3/19] Manual curation / undo-redo smoke tests
"%RDIR%\bin\Rscript.exe" "tests\unit\manual_curation_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [4/19] Stage 1 AB1 evidence helper tests
"%RDIR%\bin\Rscript.exe" "tests\unit\ab1_evidence_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [5/19] Workflow order and QC naming smoke tests
"%RDIR%\bin\Rscript.exe" "tests\contracts\workflow_order_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [6/19] Stage 2 architecture and migration tests
"%RDIR%\bin\Rscript.exe" "tests\unit\stage2_architecture_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [7/19] Stage 2 closed-gate Shiny integration tests
"%RDIR%\bin\Rscript.exe" "tests\contracts\stage2_app_contract_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [8/19] Stage 3 consensus algorithm tests
"%RDIR%\bin\Rscript.exe" "tests\unit\stage3_consensus_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [9/19] Stage 3 controlled AB1 pair tests
"%RDIR%\bin\Rscript.exe" "tests\integration\stage3_ab1_fixture_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [10/19] Stage 3 Shiny integration and startup tests
"%RDIR%\bin\Rscript.exe" "tests\integration\stage3_app_contract_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [11/19] Stage 4 multi-locus profile tests
"%RDIR%\bin\Rscript.exe" "tests\unit\stage4_multilocus_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [12/19] Stage 4 Shiny integration contract
"%RDIR%\bin\Rscript.exe" "tests\contracts\stage4_app_contract_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [13/19] Global DataTables alignment regression contract
"%RDIR%\bin\Rscript.exe" "tests\contracts\datatables_alignment_contract_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [14/19] Organized project structure contract
"%RDIR%\bin\Rscript.exe" "tests\contracts\project_structure_contract_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [15/19] Alpha 10 assay and schema 6 foundation tests
"%RDIR%\bin\Rscript.exe" "tests\unit\alpha10_assay_schema_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [16/19] Reflection follow-up correctness tests
"%RDIR%\bin\Rscript.exe" "tests\unit\reflection_followup_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [17/19] BLAST RTOE auto-poll unit tests
"%RDIR%\bin\Rscript.exe" "tests\unit\blast_polling_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [18/19] Trim/QC empty-state and amplicon overview tests
"%RDIR%\bin\Rscript.exe" "tests\unit\qc_empty_state_smoke.R"
if errorlevel 1 goto :failed

echo.
echo [19/19] Share Project snapshot tests
"%RDIR%\bin\Rscript.exe" "tests\unit\share_project_smoke.R"
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

