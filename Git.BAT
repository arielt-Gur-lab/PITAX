@echo off
setlocal
pushd "%~dp0"

powershell.exe -NoProfile -ExecutionPolicy Bypass -Command ^
"$ErrorActionPreference = 'Stop'; ^
Write-Host ''; ^
Write-Host '========================================'; ^
Write-Host '       PITAX Git Commit + Push'; ^
Write-Host '========================================'; ^
Write-Host ''; ^
git status; ^
Write-Host ''; ^
$msg = Read-Host 'Enter commit message'; ^
if ([string]::IsNullOrWhiteSpace($msg)) { ^
    Write-Host ''; ^
    Write-Host 'ERROR: Commit message is empty.' -ForegroundColor Red; ^
    Read-Host 'Press Enter to close'; ^
    exit 1 ^
}; ^
Write-Host ''; ^
Write-Host ('Commit message: ' + $msg); ^
Write-Host ''; ^
git add --all; ^
git diff --cached --quiet; ^
if ($LASTEXITCODE -eq 0) { ^
    Write-Host 'Nothing to commit.' -ForegroundColor Yellow; ^
    Read-Host 'Press Enter to close'; ^
    exit 0 ^
}; ^
Write-Host 'Committing...'; ^
git commit -m $msg; ^
if ($LASTEXITCODE -ne 0) { ^
    Write-Host ''; ^
    Write-Host 'ERROR: Commit failed.' -ForegroundColor Red; ^
    Read-Host 'Press Enter to close'; ^
    exit 1 ^
}; ^
Write-Host ''; ^
Write-Host 'Pushing...'; ^
git push; ^
if ($LASTEXITCODE -ne 0) { ^
    Write-Host ''; ^
    Write-Host 'ERROR: Push failed. The commit is saved locally.' -ForegroundColor Red; ^
    Read-Host 'Press Enter to close'; ^
    exit 1 ^
}; ^
Write-Host ''; ^
Write-Host '========================================' -ForegroundColor Green; ^
Write-Host ' COMMIT AND PUSH COMPLETED SUCCESSFULLY' -ForegroundColor Green; ^
Write-Host '========================================' -ForegroundColor Green; ^
Write-Host ''; ^
git status; ^
Write-Host ''; ^
Read-Host 'Press Enter to close'"

popd
endlocal