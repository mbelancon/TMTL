<#
  Builds the portable Windows version of TMTL2: a folder with the official
  "embeddable" Python (python.org), numpy and scipy already installed, and the program.
  Whoever receives it just unzips it and opens "Open TMTL2.bat"; no Python needed.

  Usage:  montar_portatil.bat            (at the project root)
          powershell -ExecutionPolicy Bypass -File ferramentas\montar_portatil.ps1 [-Output FOLDER] [-Python PATH]

  -Output  folder where the portable version is built (default: %USERPROFILE%\TMTL2-portable,
           outside synchronized folders such as ownCloud/OneDrive: ~150 MB and thousands of files)
  -Python  Python used only to DOWNLOAD the packages with pip (default: the user's
           Miniconda/Anaconda or the one on PATH). Any 3.x version with pip works.
#>
param(
    [Alias("Saida")][string]$Output = (Join-Path $env:USERPROFILE "TMTL2-portable"),
    [string]$Python = ""
)

$ErrorActionPreference = "Stop"
$ProgressPreference = "SilentlyContinue"          # makes Invoke-WebRequest much faster
# Python in UTF-8: with output redirected to the log, the Windows default (cp1252)
# lacks symbols such as - and chi2 and Python aborts when printing them.
$env:PYTHONUTF8 = "1"
$env:PYTHONIOENCODING = "utf-8"
try { [Console]::OutputEncoding = [Text.Encoding]::UTF8 } catch { }
[Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12

# Pinned versions: the portable version is reproducible and tested with these
$PY_VERSION = "3.12.10"
$PY_TAG     = "312"
$PACKAGES   = @("numpy==2.2.6", "scipy==1.15.3")

$root   = Split-Path -Parent $PSScriptRoot
$log    = Join-Path $root "montar_portatil.log"
Start-Transcript -Path $log -Force | Out-Null

function Step($msg) { Write-Host ""; Write-Host "== $msg" -ForegroundColor Cyan }

# Runs an external program without letting warnings on stderr (common with pip)
# stop the script; success is decided by the exit code.
function Invoke-Tool([string]$exe, [string[]]$arguments) {
    $before = $ErrorActionPreference
    $ErrorActionPreference = "Continue"
    try { & $exe @arguments 2>&1 | ForEach-Object { Write-Host $_ } }
    finally { $ErrorActionPreference = $before }
    return $LASTEXITCODE
}

try {
    $version = (Select-String -Path (Join-Path $root "pyproject.toml") -Pattern '^version\s*=\s*"(.+)"').Matches[0].Groups[1].Value
    $dest   = Join-Path $Output "TMTL2"
    $cache  = Join-Path $Output "_cache"
    $zip    = Join-Path $Output "TMTL2-$version-portable-windows.zip"
    Write-Host "TMTL2 $version - portable version"
    Write-Host "Project: $root"
    Write-Host "Output:  $Output"

    Step "1/6 Preparing folders"
    if (Test-Path $dest) { Remove-Item $dest -Recurse -Force }
    New-Item -ItemType Directory -Force -Path $dest, $cache | Out-Null

    Step "2/6 Embeddable Python $PY_VERSION (python.org)"
    $pyzip = Join-Path $cache "python-$PY_VERSION-embed-amd64.zip"
    if (-not (Test-Path $pyzip)) {
        Invoke-WebRequest "https://www.python.org/ftp/python/$PY_VERSION/python-$PY_VERSION-embed-amd64.zip" -OutFile $pyzip
    }
    $pydir = Join-Path $dest "python"
    Expand-Archive $pyzip -DestinationPath $pydir -Force
    # the ._pth file defines sys.path of the embedded Python:
    #   standard library, site-packages (numpy, scipy) and the program folder (..)
    $pth = Join-Path $pydir "python$PY_TAG._pth"
    @("python$PY_TAG.zip", ".", "Lib\site-packages", "..", "import site") | Set-Content $pth -Encoding ascii

    Step "3/6 numpy and scipy (binary packages for Windows, Python $PY_TAG)"
    if (-not $Python) {
        $candidates = @("$env:USERPROFILE\miniconda3\python.exe", "$env:USERPROFILE\anaconda3\python.exe",
                        "$env:LOCALAPPDATA\miniconda3\python.exe", "$env:LOCALAPPDATA\anaconda3\python.exe",
                        "$env:ProgramData\miniconda3\python.exe", "$env:ProgramData\anaconda3\python.exe")
        $Python = $candidates | Where-Object { Test-Path $_ } | Select-Object -First 1
        if (-not $Python) { $Python = (Get-Command python -ErrorAction SilentlyContinue).Source }
    }
    if (-not $Python) { throw "No Python with pip found to download the packages. Pass one with -Python." }
    Write-Host "Downloading with the pip of: $Python"
    $env:PATH = "$(Split-Path $Python);$(Join-Path (Split-Path $Python) 'Library\bin');$env:PATH"
    $site = Join-Path $pydir "Lib\site-packages"
    $pipArgs = @("-m", "pip", "install", "--disable-pip-version-check", "--no-warn-script-location",
                 "--target", $site, "--platform", "win_amd64", "--python-version", $PY_TAG,
                 "--implementation", "cp", "--only-binary=:all:", "--no-deps") + $PACKAGES
    if ((Invoke-Tool $Python $pipArgs) -ne 0) { throw "pip failed to download numpy/scipy" }

    Step "4/6 Copying the program"
    Copy-Item (Join-Path $root "app.py") $dest
    Copy-Item (Join-Path $root "tlfit") $dest -Recurse
    Copy-Item (Join-Path $root "tests") $dest -Recurse
    Copy-Item (Join-Path $root "README.md") $dest
    Copy-Item (Join-Path $root "README.pt-BR.md") $dest
    if (Test-Path (Join-Path $root "LICENSE")) { Copy-Item (Join-Path $root "LICENSE") $dest }
    Get-ChildItem $dest -Recurse -Directory -Filter "__pycache__" | Remove-Item -Recurse -Force

    @'
@echo off
rem TMTL2 - portable version: no installed Python needed.
rem The interface opens in the browser. To quit, close the tab (the program exits by
rem itself about 30 s later) or this window.
cd /d "%~dp0"
set PYTHONUTF8=1
"python\python.exe" app.py --browser
if errorlevel 1 pause
'@ | Set-Content (Join-Path $dest "Open TMTL2.bat") -Encoding ascii

    @'
@echo off
rem Checks that the portable TMTL2 computes correctly (takes ~30 s).
cd /d "%~dp0"
set PYTHONUTF8=1
chcp 65001 >nul
"python\python.exe" tests\validate.py
pause
'@ | Set-Content (Join-Path $dest "Check installation.bat") -Encoding ascii

    @"
TMTL2 $version - portable version for Windows
=============================================

No need to install Python or any other program.

1. Unzip this folder wherever you like (for example, in Documents).
2. Double-click "Open TMTL2.bat".
   The interface opens in the browser, in your system's language (English or
   Portuguese; there is a PT/EN switch at the top). To quit, close the tab: the
   program exits by itself about 30 s later (or close the black window).

"Check installation.bat" checks the calculations against the original code (~30 s).

Source code, documentation and new versions: https://github.com/mbelancon/TMTL
Contents: Python $PY_VERSION (python.org, embeddable), $($PACKAGES -join ', ').

License: GNU GPL version 3 or later (LICENSE file). Includes Python (PSF
license), numpy and scipy (BSD licenses), with their license texts inside the
"python" folder.

-----------------------------------------------------------------------------
Em portugues: descompacte a pasta e de dois cliques em "Open TMTL2.bat". A
interface abre no navegador, no idioma do sistema (ha um seletor PT/EN no
alto). Para encerrar, feche a aba. "Check installation.bat" confere os
calculos. Instrucoes completas em README.pt-BR.md.
"@ | Set-Content (Join-Path $dest "README.txt") -Encoding utf8

    Step "5/6 Testing the portable version"
    $pyexe = Join-Path $pydir "python.exe"
    $test = "import sys, numpy, scipy, tlfit.models, tlfit.batch; print('Python', sys.version.split()[0], '| numpy', numpy.__version__, '| scipy', scipy.__version__, '| tlfit ok')"
    Push-Location $dest
    try {
        if ((Invoke-Tool $pyexe @("-c", $test)) -ne 0) { throw "the portable version could not import the libraries" }
        if ((Invoke-Tool $pyexe @("tests\validate.py")) -ne 0) { throw "the numerical validation failed in the portable version" }
    } finally { Pop-Location }

    Step "6/6 Compressing"
    if (Test-Path $zip) { Remove-Item $zip -Force }
    Compress-Archive -Path $dest -DestinationPath $zip -CompressionLevel Optimal
    $mb = [math]::Round((Get-Item $zip).Length / 1MB, 1)
    $mbFolder = [math]::Round(((Get-ChildItem $dest -Recurse -File | Measure-Object Length -Sum).Sum) / 1MB, 1)

    Write-Host ""
    Write-Host "DONE" -ForegroundColor Green
    Write-Host "  Folder to test: $dest  ($mbFolder MB)"
    Write-Host "  Zip to distribute: $zip  ($mb MB)"
    Write-Host "  To test now: open the folder and double-click 'Open TMTL2.bat'."
    $status = 0
}
catch {
    Write-Host ""
    Write-Host "FAILED: $($_.Exception.Message)" -ForegroundColor Red
    Write-Host "The full log is in $log"
    $status = 1
}
finally {
    Stop-Transcript | Out-Null
}
exit $status
