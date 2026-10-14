param([switch]$NoBrowser)
$ErrorActionPreference = 'Stop'
$labRoot = $PSScriptRoot
$labPython = Join-Path $labRoot '.venv/Scripts/python.exe'
if (-not (Test-Path -LiteralPath $labPython)) {
    python -m venv (Join-Path $labRoot '.venv')
    if ($LASTEXITCODE -ne 0) { throw 'Could not create the Python environment.' }
    & $labPython -m pip install -r (Join-Path $labRoot 'requirements.txt')
    if ($LASTEXITCODE -ne 0) { throw 'Could not install the application dependencies.' }
}
$labAddress = 'http://127.0.0.1:8057'
$labAlreadyRunning = $false
try {
    $labResponse = Invoke-WebRequest -Uri "$labAddress/_dash-layout" -UseBasicParsing -TimeoutSec 2
    $labAlreadyRunning = $labResponse.Content.Contains('Fiber Light Lab') -or $labResponse.Content.Contains('LIGHT LAB')
} catch { }
if (-not $labAlreadyRunning) {
    $labProcess = Start-Process -FilePath $labPython -ArgumentList 'app.py' -WorkingDirectory $labRoot -WindowStyle Hidden -PassThru -RedirectStandardOutput (Join-Path $labRoot 'server.log') -RedirectStandardError (Join-Path $labRoot 'server-error.log')
    $labProcess.Id | Set-Content (Join-Path $labRoot '.server.pid')
    $labReady = $false
    for ($labAttempt = 0; $labAttempt -lt 30; $labAttempt++) {
        Start-Sleep -Milliseconds 500
        if ($labProcess.HasExited) { throw 'Application exited. See server-error.log.' }
        try { $null = Invoke-WebRequest -Uri $labAddress -UseBasicParsing -TimeoutSec 1; $labReady = $true; break } catch { }
    }
    if (-not $labReady) { throw 'Application did not become ready. See server-error.log.' }
}
if (-not $NoBrowser) { Start-Process $labAddress }
