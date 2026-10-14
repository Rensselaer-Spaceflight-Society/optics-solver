@echo off
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%~dp0LaunchDesktop.ps1"
if errorlevel 1 pause
