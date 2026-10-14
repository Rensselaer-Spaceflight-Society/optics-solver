$labPidPath = Join-Path $PSScriptRoot '.server.pid'
if (Test-Path -LiteralPath $labPidPath) {
    $labSavedId = [int](Get-Content -LiteralPath $labPidPath)
    $labProc = Get-CimInstance Win32_Process -Filter "ProcessId = $labSavedId" -ErrorAction SilentlyContinue
    $labExpected = Join-Path $PSScriptRoot '.venv\Scripts\python.exe'
    if ($labProc -and $labProc.ExecutablePath -eq $labExpected -and $labProc.CommandLine -match 'app.py') {
        Stop-Process -Id $labSavedId
        Remove-Item -LiteralPath $labPidPath
    }
}
