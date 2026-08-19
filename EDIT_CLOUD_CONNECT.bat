@echo off
cd /d "%~dp0"

REM Run the R script and wait until it finishes
"C:\Program Files\R\R-4.5.2\bin\Rscript.exe" "prepare_connect_cloud.R"

REM Open the generated manifest
start "" "manifest.json"

exit