# P-type PPR Prediction Program (PPPRPP)

A standalone web application for identifying PPR motifs in protein sequences and analyzing their binding sites on RNA targets.

## Features
- **Unified Workflow**: Input Protein and RNA sequences to perform discovery and analysis in one step.
- **Interactive Visualization**: Genome Browser with detailed binding scores and interactive site inspection.
- **Gap Analysis**: Analyze binding sites with flexible gap parameters (Gap 0, 1, 2, etc, or "Max Gap").
- **Cross-Platform**: Runs on Windows, Mac, and Linux.

## Quick Start (Mac / Linux)

### 1. Initial Setup (Run Once)

**Mac**:
- Double-click `install.command`.

**Linux**:
- Open terminal and run:
  ```bash
  bash install.command
  ```

This script will:
1.  Check/Install **HMMER** (required for motif discovery).
    -   *Mac*: Uses Homebrew (`brew install hmmer`).
    -   *Linux*: Uses apt (`sudo apt install hmmer`).
2.  Create a Python environment.
3.  Install all dependencies.

### 2. Launch App

**Mac**:
- Double-click `start.command`.

**Linux**:
- Open terminal and run:
  ```bash
  bash start.command
  ```

The application will open automatically in your default web browser.

## Troubleshooting

- **"hmmsearch not found"**: The installer attempts to install HMMER. If it fails, please install it manually:
    -   **Mac**: `brew install hmmer`
    -   **Linux**: `sudo apt install hmmer` (or within WSL on Windows)
- **Browser doesn't open**: Look for a URL in the console output (e.g., `http://127.0.0.1:7860`) and open it manually.
