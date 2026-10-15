# tasks.ps1 — convenience shortcuts for the cancer_analysis pipeline
# Usage:  .\tasks.ps1 run
#         .\tasks.ps1 test
#         .\tasks.ps1 lint
#         .\tasks.ps1 clean
#         .\tasks.ps1 help

param(
    [Parameter(Position = 0)]
    [ValidateSet("run", "test", "lint", "clean", "help")]
    [string]$Task = "help"
)

$ErrorActionPreference = "Stop"

function Invoke-Run {
    Write-Host "▶ Running pipeline..." -ForegroundColor Cyan
    python scripts/run_pipeline.py
}

function Invoke-Test {
    Write-Host "▶ Running tests..." -ForegroundColor Cyan
    pytest -q
}

function Invoke-Lint {
    Write-Host "▶ Linting with ruff..." -ForegroundColor Cyan
    ruff check src tests
}

function Invoke-Clean {
    Write-Host "▶ Cleaning build/test artifacts..." -ForegroundColor Cyan
    Remove-Item -Recurse -Force .pytest_cache -ErrorAction SilentlyContinue
    Remove-Item -Recurse -Force .ruff_cache -ErrorAction SilentlyContinue
    Get-ChildItem -Recurse -Filter __pycache__ -Directory |
        Remove-Item -Recurse -Force -ErrorAction SilentlyContinue
    Get-ChildItem -Recurse -Filter *.pyc |
        Remove-Item -Force -ErrorAction SilentlyContinue
    Get-ChildItem results -File -Exclude .gitkeep |
        Remove-Item -Force -ErrorAction SilentlyContinue
    Get-ChildItem results\figures -File -Exclude .gitkeep |
        Remove-Item -Force -ErrorAction SilentlyContinue
}

function Show-Help {
    Write-Host ""
    Write-Host "cancer-biomarker-discovery task runner" -ForegroundColor Green
    Write-Host ""
    Write-Host "  .\tasks.ps1 run     Run the full pipeline"
    Write-Host "  .\tasks.ps1 test    Run pytest"
    Write-Host "  .\tasks.ps1 lint    Run ruff"
    Write-Host "  .\tasks.ps1 clean   Remove __pycache__, test cache, and stale results"
    Write-Host "  .\tasks.ps1 help    Show this message"
    Write-Host ""
}

switch ($Task) {
    "run"  { Invoke-Run }
    "test" { Invoke-Test }
    "lint" { Invoke-Lint }
    "clean"{ Invoke-Clean }
    "help" { Show-Help }
}
