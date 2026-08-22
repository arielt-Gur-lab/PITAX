@echo off
setlocal EnableExtensions
cd /d "%~dp0"
title PITAX Git Commit + Push

echo.
echo ========================================
echo        PITAX Git Commit + Push
echo ========================================
echo.

where git.exe >nul 2>&1
if errorlevel 1 goto no_git

git.exe rev-parse --is-inside-work-tree >nul 2>&1
if errorlevel 1 goto no_repo

set "REPO_ROOT="
for /f "delims=" %%R in ('git.exe rev-parse --show-toplevel 2^>nul') do set "REPO_ROOT=%%R"
if not defined REPO_ROOT goto no_repo

cd /d "%REPO_ROOT%"
if errorlevel 1 goto no_repo

set "BRANCH="
for /f "delims=" %%B in ('git.exe branch --show-current 2^>nul') do set "BRANCH=%%B"
if /i not "%BRANCH%"=="main" goto wrong_branch

echo Repository: %REPO_ROOT%
echo Branch:     %BRANCH%
echo.

if not exist "manifest.json" goto no_manifest

set "MANIFEST_MATCHES=0"
for /f %%N in ('findstr /C:"prepare_connect_cloud.R" "manifest.json" 2^>nul ^| find /c /v ""') do set "MANIFEST_MATCHES=%%N"
if not "%MANIFEST_MATCHES%"=="1" goto stale_manifest

echo Checking origin/main before tests and commit...
git.exe fetch origin main
if errorlevel 1 goto fetch_failed

git.exe merge-base --is-ancestor origin/main HEAD >nul 2>&1
if errorlevel 1 goto remote_ahead

echo Running all PITAX tests before commit...
echo.
call "run_tests.bat"
if errorlevel 1 goto tests_failed

echo.
echo Current changes:
git.exe status --short
if errorlevel 1 goto status_failed
echo.

set "VERSION="
if exist "VERSION.txt" (
    for /f "usebackq delims=" %%V in ("VERSION.txt") do if not defined VERSION set "VERSION=%%V"
)

if defined VERSION (
    set "DEFAULT_MESSAGE=PITAX v%VERSION%"
) else (
    set "DEFAULT_MESSAGE=PITAX update"
)

set "COMMIT_MESSAGE="
set /p "COMMIT_MESSAGE=Enter commit message [%DEFAULT_MESSAGE%]: "
if not defined COMMIT_MESSAGE set "COMMIT_MESSAGE=%DEFAULT_MESSAGE%"

echo.
echo Commit message: %COMMIT_MESSAGE%
echo.

git.exe add --all
if errorlevel 1 goto add_failed

git.exe diff --cached --quiet
set "DIFF_EXIT=%ERRORLEVEL%"

if "%DIFF_EXIT%"=="0" goto push_existing
if not "%DIFF_EXIT%"=="1" goto diff_failed

echo Committing...
git.exe commit -m "%COMMIT_MESSAGE%"
if errorlevel 1 goto commit_failed

echo.
echo Pushing main to origin...
git.exe push origin main
if errorlevel 1 goto push_failed

echo.
echo ========================================
echo  COMMIT AND PUSH COMPLETED SUCCESSFULLY
echo ========================================
echo.
git.exe status

set "EXIT_CODE=0"
goto finish

:no_git
echo ERROR: Git is not installed or git.exe is not available in PATH.
goto failed

:no_repo
echo ERROR: Application is not connected to a Git repository.
echo The hidden .git folder is missing.
goto failed

:wrong_branch
echo ERROR: Expected branch main.
echo Current branch: %BRANCH%
goto failed

:no_manifest
echo ERROR: manifest.json is missing.
echo Run EDIT_CLOUD_CONNECT.bat and try again.
goto failed

:stale_manifest
echo ERROR: manifest.json still contains paths from the old structure.
echo Run EDIT_CLOUD_CONNECT.bat once and then run Git.BAT again.
goto failed

:fetch_failed
echo.
echo ERROR: Could not fetch origin/main.
echo No commit was created. Check the network or GitHub login and try again.
goto failed

:remote_ahead
echo.
echo ERROR: origin/main contains commits that are not in this local branch.
echo Nothing was committed or changed by this script.
echo Update the local repository before creating another PITAX commit.
goto failed

:tests_failed
echo.
echo ERROR: Tests failed.
echo Nothing was committed or pushed.
goto failed

:status_failed
echo ERROR: git status failed.
goto failed

:add_failed
echo ERROR: git add failed.
goto failed

:diff_failed
echo ERROR: Could not inspect the staged changes.
goto failed

:commit_failed
echo.
echo ERROR: Commit failed.
goto failed

:push_failed
echo.
echo ERROR: Push failed.
echo The commit is saved locally.
goto failed

:push_existing
echo No new file changes to commit.
echo Checking whether an existing local commit still needs to be pushed...
git.exe push origin main
if errorlevel 1 goto push_failed
echo Local and remote main are synchronized.
set "EXIT_CODE=0"
goto finish

:failed
set "EXIT_CODE=1"

:finish
echo.
pause
exit /b %EXIT_CODE%
