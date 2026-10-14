$ErrorActionPreference = 'Stop'
$labPython = Join-Path $PSScriptRoot '.venv/Scripts/python.exe'
if (-not (Test-Path -LiteralPath $labPython)) {
    python -m venv (Join-Path $PSScriptRoot '.venv')
    if ($LASTEXITCODE -ne 0) { throw 'Could not create the Python environment.' }
}
& $labPython -c 'import PySide6, matplotlib, numpy, scipy, pyqtgraph, OpenGL'
if ($LASTEXITCODE -ne 0) {
    & $labPython -m pip install -r (Join-Path $PSScriptRoot 'requirements-desktop.txt')
    if ($LASTEXITCODE -ne 0) { throw 'Could not install desktop dependencies.' }
}
& $labPython (Join-Path $PSScriptRoot 'fiber_simulation.py')
