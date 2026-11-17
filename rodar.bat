@echo off
rem Runs TMTL2 directly from the source code.
rem Looks for the user's Miniconda/Anaconda Python, then the one on PATH.
rem Optional, to open in its own window:  pip install pywebview
setlocal
cd /d "%~dp0"
set PYTHONUTF8=1
call :find_python
if not defined PY goto :no_python
rem conda folders on PATH (numpy/scipy DLLs), as the Anaconda Prompt does
for %%D in ("%PY%") do set "PATH=%%~dpD;%%~dpDLibrary\bin;%%~dpDScripts;%PATH%"
echo Using %PY%
"%PY%" app.py %*
if errorlevel 1 pause
exit /b

:find_python
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

:no_python
echo.
echo No Python with numpy and scipy was found.
echo Open the "Anaconda Prompt", go to this folder and run:  python app.py
echo (or install the libraries with:  pip install numpy scipy)
pause
exit /b 1
