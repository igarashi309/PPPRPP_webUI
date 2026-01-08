#!/bin/bash

# PPR Motif Analyzer Launcher for Mac/Linux

echo "========================================"
echo "   PPR Motif Analyzer Launcher"
echo "========================================"

# Check for Python
if command -v python3 &>/dev/null; then
    PYTHON_CMD=python3
elif command -v python &>/dev/null; then
    PYTHON_CMD=python
else
    echo "Error: Python is not installed."
    exit 1
fi

echo "Using Python: $($PYTHON_CMD --version)"

# Check for HMMER (Required for PPRFinder)
if ! command -v hmmsearch &>/dev/null; then
    echo "WARNING: 'hmmsearch' (HMMER) is not found."
    echo "PPRFinder requires HMMER."
    echo "  - Mac: brew install hmmer"
    echo "  - Linux: sudo apt install hmmer"
    read -p "Press Enter to continue anyway..."
fi

# Create Virtual Environment if not exists
VENV_DIR="venv"
if [ ! -d "$VENV_DIR" ]; then
    echo "Creating virtual environment..."
    $PYTHON_CMD -m venv $VENV_DIR
fi

# Activate Venv
source $VENV_DIR/bin/activate

# Install Dependencies
echo "Installing/Updating dependencies..."
pip install -r requirements.txt

# Setup Resources
echo "Checking resources..."
python scripts/setup_resources.py

# Run App
echo "Starting App..."
python app.py
