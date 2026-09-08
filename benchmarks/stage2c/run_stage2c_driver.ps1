param(
    [Parameter(Mandatory = $true)]
    [ValidateSet(
        "k-smoke", "methods-smoke-lambda", "methods-smoke-d",
        "methods-smoke-delta", "context-smoke", "k-formal",
        "lambda-formal", "d-formal", "delta-formal", "context-formal"
    )]
    [string]$Task,

    [string]$StagingRoot =
        "C:\Users\wenyi\AppData\Local\Temp\fastphylosig-stage2c-9825a06",

    [string]$Commit =
        "9825a0661d7a1dde6e1f0969b60d0e104291aaaf",

    [string]$OutputRoot =
        "C:\Users\wenyi\AppData\Local\Temp\fastphylosig-stage2c-9825a06\results"
)

$ErrorActionPreference = "Stop"
$rscript = "C:\Users\wenyi\AppData\Local\R\R-4.6.1\bin\Rscript.exe"
$library = Join-Path $StagingRoot "library"
$bench = Join-Path $StagingRoot "benchmarks\stage2c"
$output = Join-Path $OutputRoot $Task

if (-not (Test-Path -LiteralPath $rscript)) {
    throw "Rscript not found: $rscript"
}
if (-not (Test-Path -LiteralPath (Join-Path $library "fastphylosig\DESCRIPTION"))) {
    throw "Installed fastphylosig library not found: $library"
}

New-Item -ItemType Directory -Force -Path $output | Out-Null

$env:LC_ALL = "C"
$env:LANG = "C"
$env:LANGUAGE = "C"
$env:FASTPHYLOSIG_STAGE2C_LIBRARY = $library
$env:FASTPHYLOSIG_STAGE2C_COMMIT = $Commit
$env:FASTPHYLOSIG_STAGE2C_SOURCE_DIR = $StagingRoot
$env:FASTPHYLOSIG_STAGE2C_TIMEOUT = "300"
$env:FASTPHYLOSIG_STAGE2C_CHILD_TIMEOUT_SECONDS = "120"

switch ($Task) {
    "k-smoke" {
        & $rscript --vanilla (Join-Path $bench "run_stage2c.R") `
            $StagingRoot $output --smoke
    }
    "methods-smoke-lambda" {
        $env:FASTPHYLOSIG_STAGE2C_METHODS = "lambda"
        & $rscript --vanilla (Join-Path $bench "run_stage2c_methods.R") `
            $library $output quick
    }
    "methods-smoke-d" {
        $env:FASTPHYLOSIG_STAGE2C_METHODS = "D"
        & $rscript --vanilla (Join-Path $bench "run_stage2c_methods.R") `
            $library $output quick
    }
    "methods-smoke-delta" {
        $env:FASTPHYLOSIG_STAGE2C_METHODS = "Delta"
        & $rscript --vanilla (Join-Path $bench "run_stage2c_methods.R") `
            $library $output quick
    }
    "context-smoke" {
        & $rscript --vanilla (Join-Path $bench "run_stage2c_context.R") `
            $library $output --smoke
    }
    "k-formal" {
        & $rscript --vanilla (Join-Path $bench "run_stage2c.R") `
            $StagingRoot $output
    }
    "lambda-formal" {
        $env:FASTPHYLOSIG_STAGE2C_METHODS = "lambda"
        & $rscript --vanilla (Join-Path $bench "run_stage2c_methods.R") `
            $library $output
    }
    "d-formal" {
        $env:FASTPHYLOSIG_STAGE2C_METHODS = "D"
        & $rscript --vanilla (Join-Path $bench "run_stage2c_methods.R") `
            $library $output
    }
    "delta-formal" {
        $env:FASTPHYLOSIG_STAGE2C_METHODS = "Delta"
        & $rscript --vanilla (Join-Path $bench "run_stage2c_methods.R") `
            $library $output
    }
    "context-formal" {
        & $rscript --vanilla (Join-Path $bench "run_stage2c_context.R") `
            $library $output --formal
    }
}

if ($LASTEXITCODE -ne 0) {
    exit $LASTEXITCODE
}

Write-Output "STAGE2C_DRIVER=PASS"
Write-Output "TASK=$Task"
Write-Output "COMMIT=$Commit"
Write-Output "OUTPUT=$output"
