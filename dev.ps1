# Starts the read-only API and the Vite app. Run from the repo root.
$root = Split-Path -Parent $MyInvocation.MyCommand.Path
Write-Host "API  http://127.0.0.1:8000"
Write-Host "UI   http://127.0.0.1:5173  (mocks on unless web/.env.local sets VITE_USE_MOCKS=false)"
Start-Process -FilePath python -ArgumentList @('-m', 'uvicorn', 'api.main:app', '--host', '127.0.0.1', '--port', '8000') -WorkingDirectory $root
Start-Process -FilePath npm -ArgumentList @('run', 'dev') -WorkingDirectory (Join-Path $root 'web')
