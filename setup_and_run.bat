@echo off
setlocal

echo ========================================
echo    PPR Motif Analyzer Launcher
echo ========================================

REM Check for Python
python --version >nul 2>&1
if %errorlevel% neq 0 (
    echo Error: Python is not installed or not in PATH.
    pause
    exit /b 1
)

REM Check for HMMER (hmmsearch)
where hmmsearch >nul 2>&1
if %errorlevel% neq 0 (
    echo WARNING: 'hmmsearch' ^(HMMER^) is not found.
    echo PPRFinder requires HMMER.
    echo On Windows, you may need to install HMMER binaries and add them to PATH,
    echo or run this inside WSL.
    pause
)

REM Create Venv if missing
if not exist venv (
    echo Creating virtual environment...
    python -m venv venv
)

REM Activate Venv
call venv\Scripts\activate

REM Install Dependencies
echo Installing dependencies...
pip install -r requirements.txt

REM Setup Resources
echo Checking resources...
python scripts\setup_resources.py

REM Run App
echo Starting App...
python app.py

pause
