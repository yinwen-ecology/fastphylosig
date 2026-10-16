$ErrorActionPreference = 'Stop'
$env:LC_ALL = 'English_United States.utf8'
$env:LANG = 'English_United States.utf8'
$taskRoot = (Get-Location).Path
$taskStage = Join-Path $env:TEMP ('fpsh-' + [guid]::NewGuid().ToString('N').Substring(0,8))
$packageSource = Join-Path $taskStage 'fastphylosig'
$packageLibrary = Join-Path $taskStage 'library'
New-Item -ItemType Directory -Path $packageSource,$packageLibrary | Out-Null
foreach ($name in @('DESCRIPTION', 'NAMESPACE', 'LICENSE', 'NEWS.md', '.Rbuildignore')) {
    Copy-Item -LiteralPath (Join-Path $taskRoot $name) -Destination (Join-Path $packageSource $name)
}
foreach ($directory in @('R', 'src', 'man', 'tests', 'inst')) {
    $sourceDirectory = Join-Path $taskRoot $directory
    foreach ($file in Get-ChildItem -LiteralPath $sourceDirectory -File -Recurse) {
        if ($file.Extension -in @('.o', '.dll', '.so') -or $file.Name -eq 'symbols.rds') { continue }
        $relative = [System.IO.Path]::GetRelativePath($taskRoot, $file.FullName)
        $destination = Join-Path $packageSource $relative
        New-Item -ItemType Directory -Path (Split-Path $destination) -Force | Out-Null
        Copy-Item -LiteralPath $file.FullName -Destination $destination
    }
}
$taskStage | Set-Content -LiteralPath 'docs/signal-heterogeneity/staging-path.txt'
$packageLibrary | Set-Content -LiteralPath 'docs/signal-heterogeneity/library-path.txt'
$rExecutable = 'C:\Users\wenyi\AppData\Local\R\R-4.6.1\bin\x64\R.exe'
Push-Location $taskStage
try {
    & $rExecutable CMD build fastphylosig --no-manual *> (Join-Path $taskRoot 'docs/signal-heterogeneity/staged-build.log')
    if ($LASTEXITCODE -ne 0) { throw 'Source build failed.' }
    $archive = Join-Path $taskStage 'fastphylosig_0.3.1.9000.tar.gz'
    Copy-Item -LiteralPath $archive -Destination (Join-Path $taskRoot 'docs/signal-heterogeneity/fastphylosig_0.3.1.9000.tar.gz')
    & $rExecutable CMD INSTALL $archive "--library=$packageLibrary" --no-multiarch *> (Join-Path $taskRoot 'docs/signal-heterogeneity/install.log')
    if ($LASTEXITCODE -ne 0) { throw 'Clean source installation failed.' }
} finally {
    Pop-Location
}
