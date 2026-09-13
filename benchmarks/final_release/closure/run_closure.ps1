param([Parameter(Mandatory=$true)][string]$Stage, [switch]$ResumeInstalled)
$ErrorActionPreference = 'Stop'
$r = 'C:\Users\wenyi\AppData\Local\R\R-4.6.1\bin\R.exe'
$rs = 'C:\Users\wenyi\AppData\Local\R\R-4.6.1\bin\Rscript.exe'
$env:LC_ALL = 'C'
$env:LANG = 'C'
$env:LC_COLLATE = 'C'
$env:LC_CTYPE = 'C'
Set-Location -LiteralPath $Stage
function Run-R([string]$Exe, [string[]]$Arguments, [string]$Log) {
    $ErrorActionPreference = 'Continue'
    & $Exe @Arguments *> $Log
    $code = $LASTEXITCODE
    if ($code -ne 0) { throw "R command exited $code; see $Log" }
}
if (-not $ResumeInstalled) {
    Run-R $r @('CMD','build','source') 'build.log'
    New-Item -ItemType Directory -Path 'FINAL_AUDIT_LIB' | Out-Null
    Run-R $r @('CMD','INSTALL',('--library=' + (Join-Path $Stage 'FINAL_AUDIT_LIB')),'fastphylosig_0.2.0.tar.gz') 'install.log'
}
$env:FASTPHYLOSIG_LIBRARY = Join-Path $Stage 'FINAL_AUDIT_LIB'
$commit = (Get-Content -LiteralPath 'commit.txt' -Raw).Trim()
Run-R $rs @('--vanilla','run_scientific_validation.R','--repo',$Stage,'--out',(Join-Path $Stage 'scientific_exact_tarball'),'--source-commit',$commit) 'scientific_console.log'
Run-R $rs @('--vanilla','run_full_tests.R',(Join-Path $Stage 'testthat.Rout'),(Join-Path $Stage 'testthat_status.csv'),$commit,(Join-Path $Stage 'source')) 'full_tests_console.log'
foreach ($mode in @('no_manual','as_cran','manual')) {
    New-Item -ItemType Directory -Path $mode | Out-Null
    Push-Location $mode
    $arguments = @('CMD','check')
    if ($mode -eq 'no_manual') { $arguments += @('--no-manual','--timings') }
    if ($mode -eq 'as_cran') { $arguments += '--as-cran' }
    $arguments += (Join-Path $Stage 'fastphylosig_0.2.0.tar.gz')
    $ErrorActionPreference = 'Continue'
    & $r @arguments *> 'console.log'
    $LASTEXITCODE | Out-File -LiteralPath 'exit_code.txt'
    Pop-Location
}
