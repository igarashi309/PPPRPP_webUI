#!/bin/bash
# PPPRPP Mac Installer
# Should be runnable via Double-Click

# Ensure we are in the script directory
cd "$(dirname "$0")"

echo "=========================================="
echo "   P-type PPR Prediction Program (PPPRPP)"
echo "            Mac Installer"
echo "=========================================="
echo ""

# 1. Check/Install HMMER
echo"[1/3] Checking HMMER (hmmsearch)..."
if ! command -v hmmsearch &> /dev/null; then
    echo "   hmmsearch not found."
    
    # Detect OS
    OS="$(uname)"
    if [ "$OS" == "Darwin" ]; then
        # Mac Argument
        if command -v brew &> /dev/null; then
            echo "   Mac (Homebrew) detected. Installing hmmer..."
            brew install hmmer
        else
            echo "   Error: Homebrew (brew) not found on Mac."
            echo "   Please install it manually: brew install hmmer"
            read -p "Press Enter to exit..."
            exit 1
        fi
    elif [ "$OS" == "Linux" ]; then
        # Linux Argument
        if command -v apt-get &> /dev/null; then
            echo "   Linux (apt) detected. Attempting install (requires sudo)..."
            sudo apt-get update && sudo apt-get install -y hmmer
        else
            echo "   Error: 'apt-get' not found. Please install HMMER manually."
            read -p "Press Enter to exit..."
            exit 1
        fi
    else
        echo "   Unknown OS: $OS. Please install HMMER manually."
        read -p "Press Enter to exit..."
        exit 1
    fi
    
    # Verify installation
    if ! command -v hmmsearch &> /dev/null; then
         echo "   Error: HMMER installation failed."
         read -p "Press Enter to exit..."
         exit 1
    fi
    echo "   HMMER installed successfully."
else
    echo "   HMMER is already installed."
fi

# 2. Check Python
echo ""
echo "[2/3] Checking Python..."
if ! command -v python3 &> /dev/null; then
    echo "   Error: python3 not found."
    read -p "Press Enter to exit..."
    exit 1
fi
python3 --version

# 3. Setup Virtual Environment
echo ""
echo "[3/3] Setting up Python Environment..."

if [ ! -d "venv" ]; then
    echo "   Creating 'venv'..."
    python3 -m venv venv
else
    echo "   Using existing 'venv'."
fi

source venv/bin/activate

echo "   Installing dependencies (this may take a moment)..."
pip install --upgrade pip &> /dev/null
pip install -r requirements.txt

echo ""
echo "=========================================="
echo "        Setup Completed Successfully!"
echo "=========================================="
echo "You can now run 'start_app.command' to launch the app."
echo ""
read -p "Press Enter to close..."
