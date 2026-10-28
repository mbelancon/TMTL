<#
  Monta a versão portátil do TMTL2 para Windows: uma pasta com o Python
  "embutível" oficial (python.org), numpy e scipy já instalados e o programa.
  Quem recebe só descompacta e abre "Abrir TMTL2.bat"; não precisa ter Python.

  Uso:  montar_portatil.bat            (na raiz do projeto)
        powershell -ExecutionPolicy Bypass -File ferramentas\montar_portatil.ps1 [-Saida PASTA] [-Python CAMINHO]

  -Saida   pasta onde a versão portátil é montada (padrão: %USERPROFILE%\TMTL2-portatil,
           fora de pastas sincronizadas como ownCloud/OneDrive: são ~150 MB e milhares de arquivos)
  -Python  Python usado apenas para BAIXAR os pacotes com o pip (padrão: Miniconda/Anaconda
           do usuário ou o do PATH). Pode ser qualquer versão 3.x com pip.
#>
param(
    [string]$Saida = (Join-Path $env:USERPROFILE "TMTL2-portatil"),
    [string]$Python = ""
)

$ErrorActionPreference = "Stop"
$ProgressPreference = "SilentlyContinue"          # Invoke-WebRequest fica muito mais rápido
# Python em UTF-8: com a saída redirecionada para o registro, o padrão do Windows (cp1252)
# não tem símbolos como − e χ² e o Python aborta ao imprimi-los.
$env:PYTHONUTF8 = "1"
$env:PYTHONIOENCODING = "utf-8"
try { [Console]::OutputEncoding = [Text.Encoding]::UTF8 } catch { }
[Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12

# Versões fixas: a versão portátil é reprodutível e testada com estas
$PY_VERSION = "3.12.10"
$PY_TAG     = "312"
$PACOTES    = @("numpy==2.2.6", "scipy==1.15.3")

$raiz   = Split-Path -Parent $PSScriptRoot
$log    = Join-Path $raiz "montar_portatil.log"
Start-Transcript -Path $log -Force | Out-Null

function Etapa($msg) { Write-Host ""; Write-Host "== $msg" -ForegroundColor Cyan }

# Roda um programa externo sem que avisos no stderr (comuns no pip) interrompam o
# script; o sucesso é decidido pelo código de saída.
function Rodar([string]$exe, [string[]]$argumentos) {
    $antes = $ErrorActionPreference
    $ErrorActionPreference = "Continue"
    try { & $exe @argumentos 2>&1 | ForEach-Object { Write-Host $_ } }
    finally { $ErrorActionPreference = $antes }
    return $LASTEXITCODE
}

try {
    $versao = (Select-String -Path (Join-Path $raiz "pyproject.toml") -Pattern '^version\s*=\s*"(.+)"').Matches[0].Groups[1].Value
    $dest   = Join-Path $Saida "TMTL2"
    $cache  = Join-Path $Saida "_cache"
    $zip    = Join-Path $Saida "TMTL2-$versao-portatil-windows.zip"
    Write-Host "TMTL2 $versao - versao portatil"
    Write-Host "Projeto: $raiz"
    Write-Host "Saida:   $Saida"

    Etapa "1/6 Preparando pastas"
    if (Test-Path $dest) { Remove-Item $dest -Recurse -Force }
    New-Item -ItemType Directory -Force -Path $dest, $cache | Out-Null

    Etapa "2/6 Python $PY_VERSION embutivel (python.org)"
    $pyzip = Join-Path $cache "python-$PY_VERSION-embed-amd64.zip"
    if (-not (Test-Path $pyzip)) {
        Invoke-WebRequest "https://www.python.org/ftp/python/$PY_VERSION/python-$PY_VERSION-embed-amd64.zip" -OutFile $pyzip
    }
    $pydir = Join-Path $dest "python"
    Expand-Archive $pyzip -DestinationPath $pydir -Force
    # o arquivo ._pth define o sys.path do Python embutido:
    #   biblioteca padrão, site-packages (numpy, scipy) e a pasta do programa (..)
    $pth = Join-Path $pydir "python$PY_TAG._pth"
    @("python$PY_TAG.zip", ".", "Lib\site-packages", "..", "import site") | Set-Content $pth -Encoding ascii

    Etapa "3/6 numpy e scipy (pacotes binarios para Windows, Python $PY_TAG)"
    if (-not $Python) {
        $candidatos = @("$env:USERPROFILE\miniconda3\python.exe", "$env:USERPROFILE\anaconda3\python.exe",
                        "$env:LOCALAPPDATA\miniconda3\python.exe", "$env:LOCALAPPDATA\anaconda3\python.exe",
                        "$env:ProgramData\miniconda3\python.exe", "$env:ProgramData\anaconda3\python.exe")
        $Python = $candidatos | Where-Object { Test-Path $_ } | Select-Object -First 1
        if (-not $Python) { $Python = (Get-Command python -ErrorAction SilentlyContinue).Source }
    }
    if (-not $Python) { throw "Nao encontrei um Python com pip para baixar os pacotes. Informe com -Python." }
    Write-Host "Baixando com o pip de: $Python"
    $env:PATH = "$(Split-Path $Python);$(Join-Path (Split-Path $Python) 'Library\bin');$env:PATH"
    $site = Join-Path $pydir "Lib\site-packages"
    $pipArgs = @("-m", "pip", "install", "--disable-pip-version-check", "--no-warn-script-location",
                 "--target", $site, "--platform", "win_amd64", "--python-version", $PY_TAG,
                 "--implementation", "cp", "--only-binary=:all:", "--no-deps") + $PACOTES
    if ((Rodar $Python $pipArgs) -ne 0) { throw "o pip falhou ao baixar numpy/scipy" }

    Etapa "4/6 Copiando o programa"
    Copy-Item (Join-Path $raiz "app.py") $dest
    Copy-Item (Join-Path $raiz "tlfit") $dest -Recurse
    Copy-Item (Join-Path $raiz "tests") $dest -Recurse
    Copy-Item (Join-Path $raiz "README.md") $dest
    if (Test-Path (Join-Path $raiz "LICENSE")) { Copy-Item (Join-Path $raiz "LICENSE") $dest }
    Get-ChildItem $dest -Recurse -Directory -Filter "__pycache__" | Remove-Item -Recurse -Force

    @'
@echo off
rem TMTL2 - versao portatil: nao precisa de Python instalado.
rem A interface abre no navegador. Para encerrar, feche a aba (o programa sai sozinho
rem cerca de 30 s depois) ou esta janela.
cd /d "%~dp0"
set PYTHONUTF8=1
"python\python.exe" app.py --browser
if errorlevel 1 pause
'@ | Set-Content (Join-Path $dest "Abrir TMTL2.bat") -Encoding ascii

    @'
@echo off
rem Confere se o TMTL2 portatil calcula corretamente (leva ~30 s).
cd /d "%~dp0"
set PYTHONUTF8=1
chcp 65001 >nul
"python\python.exe" tests\validate.py
pause
'@ | Set-Content (Join-Path $dest "Verificar instalacao.bat") -Encoding ascii

    @"
TMTL2 $versao - versao portatil para Windows
============================================

Nao precisa instalar Python nem nenhum outro programa.

1. Descompacte esta pasta onde quiser (por exemplo, em Documentos).
2. De dois cliques em "Abrir TMTL2.bat".
   A interface abre no navegador. Para encerrar, feche a aba: o programa
   sai sozinho cerca de 30 s depois (ou feche a janela preta).

"Verificar instalacao.bat" confere os calculos contra o codigo original (~30 s).

Codigo-fonte, documentacao e novas versoes: https://github.com/mbelancon/TMTL
Conteudo: Python $PY_VERSION (python.org, embutivel), $($PACOTES -join ', ').

Licenca: GNU GPL versao 3 ou posterior (arquivo LICENSE). Inclui o Python
(licenca PSF), o numpy e o scipy (licencas BSD), com seus textos de licenca
dentro da pasta "python".
"@ | Set-Content (Join-Path $dest "LEIA-ME.txt") -Encoding utf8

    Etapa "5/6 Testando a versao portatil"
    $pyexe = Join-Path $pydir "python.exe"
    $teste = "import sys, numpy, scipy, tlfit.models, tlfit.batch; print('Python', sys.version.split()[0], '| numpy', numpy.__version__, '| scipy', scipy.__version__, '| tlfit ok')"
    Push-Location $dest
    try {
        if ((Rodar $pyexe @("-c", $teste)) -ne 0) { throw "a versao portatil nao conseguiu importar as bibliotecas" }
        if ((Rodar $pyexe @("tests\validate.py")) -ne 0) { throw "a validacao numerica falhou na versao portatil" }
    } finally { Pop-Location }

    Etapa "6/6 Compactando"
    if (Test-Path $zip) { Remove-Item $zip -Force }
    Compress-Archive -Path $dest -DestinationPath $zip -CompressionLevel Optimal
    $mb = [math]::Round((Get-Item $zip).Length / 1MB, 1)
    $mbPasta = [math]::Round(((Get-ChildItem $dest -Recurse -File | Measure-Object Length -Sum).Sum) / 1MB, 1)

    Write-Host ""
    Write-Host "PRONTO" -ForegroundColor Green
    Write-Host "  Pasta para testar: $dest  ($mbPasta MB)"
    Write-Host "  Zip para distribuir: $zip  ($mb MB)"
    Write-Host "  Para testar agora: abra a pasta e de dois cliques em 'Abrir TMTL2.bat'."
    $status = 0
}
catch {
    Write-Host ""
    Write-Host "FALHOU: $($_.Exception.Message)" -ForegroundColor Red
    Write-Host "O registro completo esta em $log"
    $status = 1
}
finally {
    Stop-Transcript | Out-Null
}
exit $status
