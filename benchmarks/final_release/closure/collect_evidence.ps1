param([Parameter(Mandatory=$true)][string]$Stage,
      [Parameter(Mandatory=$true)][string]$Repository)
$ErrorActionPreference = 'Stop'
$destination = Join-Path $Repository 'benchmarks/final_release/closure/results'
New-Item -ItemType Directory -Force -Path $destination | Out-Null
foreach ($name in @('build.log','install.log','scientific_console.log',
                    'testthat.Rout','testthat_status.csv','full_tests_console.log',
                    'commit.txt','scientific_harness_path_failure_console.log')) {
    Copy-Item -LiteralPath (Join-Path $Stage $name) -Destination $destination
}
foreach ($name in @('scientific_exact_tarball','scientific_harness_path_failure',
                    'no_manual','as_cran','manual','extra_gates')) {
    $source = Join-Path $Stage $name
    if (Test-Path -LiteralPath $source) {
        $target = Join-Path $destination $name
        New-Item -ItemType Directory -Force -Path $target | Out-Null
        if ($name -in @('no_manual','as_cran','manual')) {
            Copy-Item -LiteralPath (Join-Path $source 'console.log'),(Join-Path $source 'exit_code.txt') -Destination $target
            Copy-Item -LiteralPath (Join-Path $source 'fastphylosig.Rcheck/00check.log') -Destination $target
        } else {
            Copy-Item -Path (Join-Path $source '*') -Destination $target -Recurse
        }
    }
}
Copy-Item -LiteralPath (Join-Path $Stage 'fastphylosig_0.2.0.tar.gz') -Destination $Repository -Force
$entries = Get-ChildItem -LiteralPath $destination -File -Recurse | ForEach-Object {
    [pscustomobject]@{path=$_.FullName.Substring($Repository.Length + 1).Replace('\','/'); sha256=(Get-FileHash -LiteralPath $_.FullName -Algorithm SHA256).Hash.ToLower(); bytes=$_.Length}
}
$entries | Export-Csv -LiteralPath (Join-Path $destination 'evidence_sha256.csv') -NoTypeInformation -Encoding UTF8
