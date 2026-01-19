# P-type PPR Prediction Program (PPPRPP)

A standalone web application for identifying PPR motifs in protein sequences and analyzing their binding sites on RNA targets.

## Features
- **Unified Workflow**: Input Protein and RNA sequences to perform discovery and analysis in one step.
- **Interactive Visualization**: Genome Browser with detailed binding scores and interactive site inspection.
- **Gap Analysis**: Analyze binding sites with specific Gap Sizes (e.g., Gap 3). The analysis efficiently scans for sites matching the specified gap size (in addition to Gap 0). When inspecting a specific site, you can interactively explore all variations (e.g., Gap 1, 2, 3) instantly.
- **Cross-Platform**: Runs on Windows (via WSL), Mac, and Linux.

## Quick Start (Mac / Linux / Windows WSL)

### 1. Initial Setup (Run Once)

**Mac**:
- Double-click `install.command`.

**Linux / Windows (WSL)**:
- **Windows Users**: Please enable **WSL (Windows Subsystem for Linux)** and install Ubuntu. Run all commands inside the WSL terminal.
- Open terminal and run:
  ```bash
  bash install.command
  ```

This script will:
1.  Check/Install **HMMER** (required for motif discovery).
    -   *Mac*: Uses Homebrew (`brew install hmmer`).
    -   *Linux/WSL*: Uses apt (`sudo apt install hmmer`).
2.  Create a Python environment.
3.  Install all dependencies.

### 2. Launch App

**Mac**:
- Double-click `start.command`.

**Linux / Windows (WSL)**:
- Open terminal and run:
  ```bash
  bash start.command
  ```

The application will open automatically in your default web browser.

## Troubleshooting

- **"hmmsearch not found"**: The installer attempts to install HMMER. If it fails, please install it manually:
    -   **Mac**: `brew install hmmer`
    -   **Linux/WSL**: `sudo apt install hmmer`
- **Browser doesn't open**: Look for a URL in the console output (e.g., `http://127.0.0.1:7860`) and open it manually.
