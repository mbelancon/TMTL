@echo off
rem Roda o TMTL2 direto do codigo (sem gerar executavel).
rem Procura o Python do Miniconda/Anaconda do usuario, depois o do PATH.
rem Opcional, para abrir em janela propria:  pip install pywebview
setlocal
cd /d "%~dp0"
set PYTHONUTF8=1
call :achar_python
if not defined PY goto :sem_python
rem pastas do conda no PATH (DLLs do numpy/scipy), como faz o Anaconda Prompt
for %%D in ("%PY%") do set "PATH=%%~dpD;%%~dpDLibrary\bin;%%~dpDScripts;%PATH%"
echo Usando %PY%
"%PY%" app.py %*
if errorlevel 1 pause
exit /b

:achar_python
for %%P in (
    "%USERPROFILE%\miniconda3\python.exe"
    "%USERPROFILE%\anaconda3\python.exe"
    "%LOCALAPPDATA%\miniconda3\python.exe"
    "%LOCALAPPDATA%\anaconda3\python.exe"
    "%ProgramData%\miniconda3\python.exe"
    "%ProgramData%\anaconda3\python.exe"
) do (
    if not defined PY if exist "%%~P" set "PY=%%~P"
)
if not defined PY (
    python -c "import numpy, scipy" >nul 2>&1 && set "PY=python"
)
if not defined PY (
    py -3 -c "import numpy, scipy" >nul 2>&1 && set "PY=py"
)
exit /b

:sem_python
echo.
echo Nao encontrei um Python com numpy e scipy.
echo Abra o "Anaconda Prompt", va ate esta pasta e rode:  python app.py
echo (ou instale as bibliotecas com:  pip install numpy scipy)
pause
exit /b 1
