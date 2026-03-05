#!/bin/bash
# PPPRPP Mac Launcher
# Should be runnable via Double-Click

# Ensure we are in the script directory
cd "$(dirname "$0")"

# Check if venv exists
if [ ! -d "venv" ]; then
    echo "Error: 'venv' directory not found."
    echo "Please run 'install.command' first to set up the environment."
    echo ""
    read -p "Press Enter to exit..."
    exit 1
fi

# Activate environment and run
source venv/bin/activate

echo "Starting PPPRPP..."
echo "The application should open in your browser shortly."
echo "Keep this window open while using the app."
echo ""

python3 app.py

# If app exits (e.g. crash), pause so user can read error
echo ""
echo "Application stopped."
read -p "Press Enter to close..."
