<#
.SYNOPSIS
    Reorganiza X:\AnaGio em um repositorio publicavel, separando codigo de dados.

.DESCRIPTION
    Por padrao roda em SIMULACAO: mostra o que faria e nao altera nada.
    Para aplicar de verdade, use -Executar.

    Principios:
      - Dado de paciente nunca e apagado, apenas MOVIDO para fora do repositorio.
      - A .venv NAO e removida: ela e necessaria para trabalhar e o .gitignore
        ja impede que suba para o git. Apagar so daria retrabalho.
      - Duplicatas na raiz sao movidas para src/kirsolve/, sobrescrevendo a
        versao antiga (a da raiz e a mais nova).

.EXAMPLE
    .\organizar_repo.ps1
    Mostra o plano sem alterar nada.

.EXAMPLE
    .\organizar_repo.ps1 -Executar
    Aplica as mudancas.

.EXAMPLE
    .\organizar_repo.ps1 -Raiz "X:\AnaGio" -PastaDados "X:\AnaGio_dados" -Executar
#>

[CmdletBinding()]
param(
    [string]$Raiz = "X:\AnaGio",
    [string]$PastaDados = "",
    [switch]$Executar
)

$ErrorActionPreference = "Stop"

if (-not (Test-Path $Raiz)) {
    Write-Host "ERRO: nao encontrei $Raiz" -ForegroundColor Red
    exit 1
}
$Raiz = (Resolve-Path $Raiz).Path
if ([string]::IsNullOrWhiteSpace($PastaDados)) {
    $PastaDados = Join-Path (Split-Path $Raiz -Parent) ((Split-Path $Raiz -Leaf) + "_dados")
}

$modo = if ($Executar) { "EXECUTANDO" } else { "SIMULACAO (use -Executar para aplicar)" }
Write-Host ""
Write-Host "=== Organizacao do repositorio kirsolve ===" -ForegroundColor Cyan
Write-Host "Repositorio : $Raiz"
Write-Host "Dados vao p/: $PastaDados"
Write-Host "Modo        : $modo" -ForegroundColor Yellow
Write-Host ""

$acoes = @()

function Add-Acao($tipo, $de, $para, $motivo) {
    $script:acoes += [pscustomobject]@{
        Tipo = $tipo; De = $de; Para = $para; Motivo = $motivo
    }
}

# ---------------------------------------------------------------- 1. duplicatas
$duplicatas = @("cli.py", "loaders.py", "validate.py")
foreach ($f in $duplicatas) {
    $origem = Join-Path $Raiz $f
    if (Test-Path $origem) {
        $destino = Join-Path $Raiz "src\kirsolve\$f"
        $nota = "versao nova solta na raiz; o pacote instalado le de src/"
        if (Test-Path $destino) {
            $tamOrigem  = (Get-Item $origem).Length
            $tamDestino = (Get-Item $destino).Length
            $nota = "$nota (raiz $tamOrigem B -> src $tamDestino B)"
        }
        Add-Acao "MOVER" $f "src\kirsolve\$f" $nota
    }
}

# ---------------------------------------------------------------- 2. dados
$padroesDados = @("*.xlsx", "*.xls", "*.bam", "*.cram", "*.fq.gz", "*.fastq.gz")
foreach ($padrao in $padroesDados) {
    Get-ChildItem -Path $Raiz -Filter $padrao -File -ErrorAction SilentlyContinue | ForEach-Object {
        Add-Acao "MOVER P/ FORA" $_.Name "$PastaDados\$($_.Name)" "dado ou resultado; nao vai para repositorio publico"
    }
}
foreach ($pasta in @("resultados", "resultados2", "resultados_kirsolve", "benchmark", "dados", "data")) {
    $p = Join-Path $Raiz $pasta
    if (Test-Path $p) {
        Add-Acao "MOVER P/ FORA" "$pasta\" "$PastaDados\$pasta\" "saidas geradas a partir de dados reais"
    }
}

# ---------------------------------------------------------------- 3. descartaveis
$descartaveis = @{
    "benchmark_hprc.sh" = "substituido pelo RUNBOOK; continha sintaxe incorreta das ferramentas"
}
foreach ($f in $descartaveis.Keys) {
    if (Test-Path (Join-Path $Raiz $f)) {
        Add-Acao "APAGAR" $f "-" $descartaveis[$f]
    }
}

# pastas vazias
foreach ($pasta in @("config")) {
    $p = Join-Path $Raiz $pasta
    if (Test-Path $p) {
        $n = @(Get-ChildItem $p -Force -ErrorAction SilentlyContinue).Count
        if ($n -eq 0) { Add-Acao "APAGAR" "$pasta\" "-" "pasta vazia" }
    }
}

# ---------------------------------------------------------------- 3b. scripts soltos
$scriptsSoltos = @("immuannot_para_verdade.py", "gerar_exemplo.py", "slurm_kirsolve.sh")
foreach ($f in $scriptsSoltos) {
    if (Test-Path (Join-Path $Raiz $f)) {
        Add-Acao "MOVER" $f "scripts\$f" "script auxiliar; o lugar dele e scripts/"
    }
}

# ---------------------------------------------------------------- 3c. lixo de build
foreach ($p in @(Get-ChildItem -Path $Raiz -Recurse -Directory -Force -ErrorAction SilentlyContinue |
                 Where-Object { $_.FullName -notmatch '\\\.venv\\' -and
                                ($_.Name -eq "__pycache__" -or $_.Name -like "*.egg-info" -or
                                 $_.Name -eq ".pytest_cache") })) {
    $rel = $p.FullName.Replace("$Raiz\", "")
    Add-Acao "APAGAR" "$rel\" "-" "gerado automaticamente; o pip recria"
}

# ---------------------------------------------------------------- 4. destinos fixos
# Arquivos baixados um a um caem na raiz. Sem isso, o README referencia
# assets/logo.jpeg e .github/workflows/ci.yml que nunca existiram, e os
# links quebram no GitHub.
$destinos = @{
    "logo.jpeg"                 = "assets"
    "Gemini_Generated_Image_5ncme85ncme85ncm.jpeg" = "assets"
    "ci.yml"                    = ".github\workflows"
    "pages.yml"                 = ".github\workflows"
    "index.html"                = "docs"
    "validacao.md"              = "docs"
    "cloud-config-hetzner.yaml" = "docs"
    "RUNBOOK_servidor.md"       = "docs"
    "launch.json"               = ".vscode"
    "settings.json"             = ".vscode"
    "test_kirsolve.py"          = "tests"
    "nomenclature.py"           = "src\kirsolve"
    "parsing.py"                = "src\kirsolve"
    "loaders.py"                = "src\kirsolve"
    "em.py"                     = "src\kirsolve"
    "resolve.py"                = "src\kirsolve"
    "priors.py"                 = "src\kirsolve"
    "report.py"                 = "src\kirsolve"
    "validate.py"               = "src\kirsolve"
    "cli.py"                    = "src\kirsolve"
}
foreach ($f in $destinos.Keys) {
    if (Test-Path (Join-Path $Raiz $f)) {
        $pasta = $destinos[$f]
        $nomeFinal = if ($f -like "Gemini_Generated*") { "logo.jpeg" } else { $f }
        Add-Acao "MOVER" $f "$pasta\$nomeFinal" "o README e o codigo esperam esse caminho"
    }
}

# ---------------------------------------------------------------- plano
if ($acoes.Count -eq 0) {
    Write-Host "Nada a fazer: a pasta ja esta organizada." -ForegroundColor Green
} else {
    $acoes | Format-Table Tipo, De, Para, Motivo -AutoSize -Wrap
}

# ---------------------------------------------------------------- executar
if ($Executar -and $acoes.Count -gt 0) {
    Write-Host "Aplicando..." -ForegroundColor Yellow

    New-Item -ItemType Directory -Force -Path $PastaDados | Out-Null
    New-Item -ItemType Directory -Force -Path (Join-Path $Raiz "src\kirsolve") | Out-Null
    New-Item -ItemType Directory -Force -Path (Join-Path $Raiz "docs") | Out-Null

    foreach ($a in $acoes) {
        $origem = Join-Path $Raiz ($a.De -replace '\\$', '')
        try {
            switch ($a.Tipo) {
                "MOVER" {
                    $destino = Join-Path $Raiz ($a.Para -replace '\\$', '')
                    $pai = Split-Path $destino -Parent
                    New-Item -ItemType Directory -Force -Path $pai | Out-Null
                    Move-Item -Path $origem -Destination $destino -Force
                    Write-Host "  movido   $($a.De)" -ForegroundColor Green
                }
                "MOVER P/ FORA" {
                    $destino = $a.Para -replace '\\$', ''
                    $pai = Split-Path $destino -Parent
                    New-Item -ItemType Directory -Force -Path $pai | Out-Null
                    Move-Item -Path $origem -Destination $destino -Force
                    Write-Host "  para fora $($a.De)" -ForegroundColor Cyan
                }
                "APAGAR" {
                    Remove-Item -Path $origem -Recurse -Force
                    Write-Host "  apagado  $($a.De)" -ForegroundColor DarkGray
                }
            }
        } catch {
            Write-Host "  FALHOU   $($a.De): $_" -ForegroundColor Red
        }
    }
    Write-Host ""
}

# ---------------------------------------------------------------- conferencia
Write-Host "=== Conferencia ===" -ForegroundColor Cyan

$obrigatorios = @(
    "pyproject.toml", "README.md", "CHANGELOG.md", "LICENSE", "CITATION.cff",
    ".gitignore", "Dockerfile", "Makefile", "app.py",
    "assets\logo.jpeg", "docs\index.html", "docs\validacao.md",
    ".github\workflows\ci.yml", ".github\workflows\pages.yml",
    "src\kirsolve\__init__.py", "src\kirsolve\cli.py", "src\kirsolve\loaders.py",
    "src\kirsolve\validate.py", "src\kirsolve\resolve.py", "src\kirsolve\em.py",
    "src\kirsolve\nomenclature.py", "src\kirsolve\parsing.py",
    "src\kirsolve\priors.py", "src\kirsolve\report.py",
    "tests\test_kirsolve.py",
    "scripts\gerar_exemplo.py", "scripts\immuannot_para_verdade.py"
)
$faltando = @()
foreach ($f in $obrigatorios) {
    if (-not (Test-Path (Join-Path $Raiz $f))) { $faltando += $f }
}
if ($faltando.Count -eq 0) {
    Write-Host "  Todos os arquivos do repositorio estao presentes." -ForegroundColor Green
} else {
    Write-Host "  Faltando (baixe da conversa e coloque no lugar):" -ForegroundColor Yellow
    $faltando | ForEach-Object { Write-Host "    $_" -ForegroundColor Yellow }
}

# links do README apontam para arquivos que existem?
Write-Host ""
Write-Host "  Conferindo os links do README..." -NoNewline
$readme = Join-Path $Raiz "README.md"
$linksQuebrados = @()
if (Test-Path $readme) {
    $texto = Get-Content $readme -Raw
    foreach ($m in [regex]::Matches($texto, '\]\(([^)]+)\)')) {
        $alvo = $m.Groups[1].Value
        if ($alvo -match '^(https?://|#)') { continue }
        if (-not (Test-Path (Join-Path $Raiz ($alvo -replace '/', '\')))) {
            $linksQuebrados += $alvo
        }
    }
    foreach ($m in [regex]::Matches($texto, '<img src="([^"]+)"')) {
        $alvo = $m.Groups[1].Value
        if ($alvo -match '^https?://') { continue }
        if (-not (Test-Path (Join-Path $Raiz ($alvo -replace '/', '\')))) {
            $linksQuebrados += "imagem: $alvo"
        }
    }
}
if ($linksQuebrados.Count -eq 0) {
    Write-Host " todos ok." -ForegroundColor Green
} else {
    Write-Host ""
    Write-Host "  LINKS QUEBRADOS no README (vao dar 404 no GitHub):" -ForegroundColor Red
    $linksQuebrados | Sort-Object -Unique | ForEach-Object { Write-Host "    $_" -ForegroundColor Red }
}

# varredura por dado sensivel remanescente
Write-Host ""
Write-Host "  Varredura por dado identificavel..." -NoNewline
$suspeitos = @()
Get-ChildItem -Path $Raiz -Recurse -File -Include *.py,*.md,*.toml,*.yaml,*.yml,*.cff,*.txt `
    -ErrorAction SilentlyContinue |
    Where-Object { $_.FullName -notmatch '\\\.venv\\' } |
    ForEach-Object {
        $conteudo = Get-Content $_.FullName -Raw -ErrorAction SilentlyContinue
        if ($conteudo -and ($conteudo -match 'ART00\d|Giovana|49\.12\.97\.7')) {
            $suspeitos += $_.FullName.Replace("$Raiz\", "")
        }
    }
if ($suspeitos.Count -eq 0) {
    Write-Host " limpo." -ForegroundColor Green
} else {
    Write-Host ""
    Write-Host "  ATENCAO - identificadores encontrados em:" -ForegroundColor Red
    $suspeitos | ForEach-Object { Write-Host "    $_" -ForegroundColor Red }
}

# ---------------------------------------------------------------- proximos passos
Write-Host ""
if (-not $Executar) {
    Write-Host "Nada foi alterado. Para aplicar:" -ForegroundColor Yellow
    Write-Host "  .\organizar_repo.ps1 -Executar"
} else {
    Write-Host "Pronto. Proximos passos:" -ForegroundColor Green
    Write-Host "  cd `"$Raiz`""
    Write-Host "  pip install -e ."
    Write-Host "  pytest -q                          # esperado: 20 passed"
    Write-Host "  python scripts\gerar_exemplo.py exemplo --amostras 60"
    Write-Host "  python -m kirsolve benchmark exemplo\ping exemplo\kirmapper --truth exemplo\verdade.tsv -o benchmark --sample-regex '^(SIM\d+)'"
    Write-Host ""
    Write-Host "  git init; git add .; git status    # confira ANTES de commitar"
    Write-Host ""
    Write-Host "  A .venv continua no lugar: e necessaria para trabalhar e o" -ForegroundColor DarkGray
    Write-Host "  .gitignore ja impede que suba. Seus dados estao em:" -ForegroundColor DarkGray
    Write-Host "  $PastaDados" -ForegroundColor DarkGray
}
Write-Host ""
