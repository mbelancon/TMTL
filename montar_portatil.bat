@echo off
rem Monta a versao portatil do TMTL2 (Python embutido + numpy + scipy + programa).
rem Precisa de internet. O resultado fica em %USERPROFILE%\TMTL2-portatil
rem (fora do ownCloud/OneDrive: sao ~150 MB e milhares de arquivos).
cd /d "%~dp0"
powershell -NoProfile -ExecutionPolicy Bypass -File "ferramentas\montar_portatil.ps1" %*
pause
