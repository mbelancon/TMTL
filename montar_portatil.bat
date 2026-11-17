@echo off
rem Builds the portable version of TMTL2 (embedded Python + numpy + scipy + program).
rem Needs internet. The result goes to %USERPROFILE%\TMTL2-portable
rem (outside ownCloud/OneDrive: ~150 MB and thousands of files).
cd /d "%~dp0"
powershell -NoProfile -ExecutionPolicy Bypass -File "ferramentas\montar_portatil.ps1" %*
pause
