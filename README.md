# PPR Motif Analyzer (Local Version)

A standalone web application for identifying PPR motifs in protein sequences and analyzing their binding sites on RNA targets.

## Features
- **Unified Workflow**: Input Protein and RNA sequences to perform discovery and analysis in one step.
- **Interactive Visualization**: Genome Browser with clickable sites to view detailed binding scores.
- **Cross-Platform**: Runs on Windows, Mac, and Linux.

## Prerequisites

Before running, ensure you have the following installed:

1.  **Python 3.10+**
2.  **HMMER (`hmmsearch`)**: Required for motif discovery.
    - **Mac**: `brew install hmmer`
    - **Linux**: `sudo apt install hmmer`
    - **Windows**: Download HMMER binaries (http://hmmer.org/) and add to your system PATH, or run this tool inside WSL (Windows Subsystem for Linux).
    - **Mac**: `brew install hmmer`
    - **Linux**: `sudo apt install hmmer`
    - **Windows**: Download HMMER binaries (http://hmmer.org/) and add to your system PATH, or run this tool inside WSL (Windows Subsystem for Linux).

## How to Run

### Mac / Linux
1.  Open Terminal.
2.  Navigate to this folder.
3.  Run the launch script:
    ```bash
    bash setup_and_run.sh
    ```

### Windows
1.  Double-click `setup_and_run.bat`.

The script will automatically:
1.  Create a Python virtual environment (`venv`).
2.  Install necessary libraries (`gradio`, `biopython`, etc.).
3.  Download required data files (`PPRfinder.jar`, HMM models).
4.  Launch the application in your default web browser.

## Troubleshooting

- **"Java not found"**: Install Java and ensure `java` command works in terminal/cmd.
- **"hmmsearch not found"**: Install HMMER. On Windows, this is the most common issue; using WSL is recommended if native binaries are hard to configure.
- **Browser doesn't open**: Look for a URL in the console output (e.g., `http://127.0.0.1:7860`) and open it manually.
