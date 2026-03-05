#!/bin/bash
# PPPRPP Mac Installer
# Should be runnable via Double-Click

# Ensure we are in the script directory
cd "$(dirname "$0")"
OS="$(uname)"

get_java_major() {
    local java_bin="$1"
    local raw
    raw="$("$java_bin" -version 2>&1 | awk -F '\"' '/version/ {print $2}' | head -n 1)"
    if [[ -z "$raw" ]]; then
        echo ""
        return
    fi
    if [[ "$raw" == 1.* ]]; then
        echo "$raw" | awk -F '.' '{print $2}'
    else
        echo "$raw" | awk -F '.' '{print $1}'
    fi
}

echo "=========================================="
echo "   P-type PPR Prediction Program (PPPRPP)"
echo "            Mac Installer"
echo "=========================================="
echo ""

# 1. Check/Install HMMER
echo "[1/4] Checking HMMER (hmmsearch)..."
if ! command -v hmmsearch &> /dev/null; then
    echo "   hmmsearch not found."
    
    # Detect OS
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

# 2. Check/Install Java 17+
echo ""
echo "[2/4] Checking Java (17+) for PPRFinder..."
JAVA_OK=0
JAVA_BIN=""

if command -v java &> /dev/null; then
    JAVA_MAJOR="$(get_java_major "$(command -v java)")"
    if [[ -n "$JAVA_MAJOR" && "$JAVA_MAJOR" -ge 17 ]]; then
        JAVA_OK=1
        JAVA_BIN="$(command -v java)"
    fi
fi

if [ "$JAVA_OK" -eq 0 ]; then
    for CAND in "/opt/homebrew/opt/openjdk@17/bin/java" "/usr/local/opt/openjdk@17/bin/java" "/usr/lib/jvm/java-17-openjdk-amd64/bin/java" "/usr/lib/jvm/java-17-openjdk/bin/java"; do
        if [ -x "$CAND" ]; then
            JAVA_OK=1
            JAVA_BIN="$CAND"
            break
        fi
    done
fi

if [ "$JAVA_OK" -eq 0 ]; then
    if [ "$OS" == "Darwin" ]; then
        if command -v brew &> /dev/null; then
            echo "   Java 17 not found. Installing openjdk@17..."
            brew install openjdk@17
        else
            echo "   Error: Homebrew (brew) not found on Mac."
            echo "   Please install Java 17 manually."
            read -p "Press Enter to exit..."
            exit 1
        fi
    elif [ "$OS" == "Linux" ]; then
        if command -v apt-get &> /dev/null; then
            echo "   Java 17 not found. Attempting install (requires sudo)..."
            sudo apt-get update && sudo apt-get install -y openjdk-17-jre-headless
        else
            echo "   Error: 'apt-get' not found. Please install Java 17 manually."
            read -p "Press Enter to exit..."
            exit 1
        fi
    else
        echo "   Unknown OS: $OS. Please install Java 17 manually."
        read -p "Press Enter to exit..."
        exit 1
    fi

    if command -v java &> /dev/null; then
        JAVA_MAJOR="$(get_java_major "$(command -v java)")"
        if [[ -n "$JAVA_MAJOR" && "$JAVA_MAJOR" -ge 17 ]]; then
            JAVA_OK=1
            JAVA_BIN="$(command -v java)"
        fi
    fi
    if [ "$JAVA_OK" -eq 0 ]; then
        for CAND in "/opt/homebrew/opt/openjdk@17/bin/java" "/usr/local/opt/openjdk@17/bin/java" "/usr/lib/jvm/java-17-openjdk-amd64/bin/java" "/usr/lib/jvm/java-17-openjdk/bin/java"; do
            if [ -x "$CAND" ]; then
                JAVA_OK=1
                JAVA_BIN="$CAND"
                break
            fi
        done
    fi
fi

if [ "$JAVA_OK" -eq 1 ]; then
    echo "   Java 17+ found: $JAVA_BIN"
else
    echo "   Error: Java 17+ not found."
    echo "   Please install Java 17 and retry."
    read -p "Press Enter to exit..."
    exit 1
fi

# 3. Check Python
echo ""
echo "[3/4] Checking Python..."
if ! command -v python3 &> /dev/null; then
    echo "   Error: python3 not found."
    echo "   Please install Python 3.10 or newer from https://www.python.org/downloads/"
    read -p "Press Enter to exit..."
    exit 1
fi
python3 --version
PY_MINOR=$(python3 -c "import sys; print(sys.version_info.minor)")
PY_MAJOR=$(python3 -c "import sys; print(sys.version_info.major)")
if [[ "$PY_MAJOR" -lt 3 || ( "$PY_MAJOR" -eq 3 && "$PY_MINOR" -lt 10 ) ]]; then
    echo "   Error: Python 3.10 or newer is required (found $(python3 --version))."
    echo "   gradio 6.x requires Python >=3.10."
    echo "   Please install Python 3.10+ from https://www.python.org/downloads/"
    read -p "Press Enter to exit..."
    exit 1
fi
echo "   Python version OK (>=3.10)."

# 4. Setup Virtual Environment
echo ""
echo "[4/4] Setting up Python Environment..."

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
echo "You can now run 'start.command' to launch the app."
echo ""
read -p "Press Enter to close..."
