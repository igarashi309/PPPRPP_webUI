# PPR Target Site Prediction Program (PPPRPP)

A standalone web application for integrated PPR motif discovery and RNA target site prediction. The tool accepts a PPR protein sequence and one or more RNA sequences as input, automatically identifies PPR motifs, scores all candidate binding windows, and visualizes the results—no bioinformatics expertise required. Both P-type and PLS-type scoring schemes are supported.

For motif discovery, the app uses **PPRFinder** (Java implementation, `PPRfinder.jar`) or the **ScanProsite API** (network-dependent; see notes below).

## Features
- **Integrated Workflow**: Input protein and RNA sequences to perform motif discovery, target scoring, and visualization in one step.
- **P-type and PLS-type Support**: Select the scoring method appropriate for your PPR protein.
- **Interactive Visualization**: Genome Browser with binding score distribution and interactive site inspection.
- **Gap Analysis**: Analyze binding sites with a specified Gap Size (e.g., Gap 3). The analysis evaluates the specified gap size together with Gap 0 as a baseline. When inspecting a specific site, you can interactively explore all gap variations (e.g., Gap 1, 2, 3) instantly.
- **Strand Filtering**: The **Show Strand** option filters the browser display to Forward-only, Reverse-only, or Both strands. When "Both" is selected, forward and reverse scores appear in separate subplots. The current strand setting is shown in the plot title.
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

**Requirements**: Python 3.10 or newer is required (gradio 6.x dependency).

This script will:
1.  Check/Install **HMMER** (required for motif discovery).
    -   *Mac*: Uses Homebrew (`brew install hmmer`).
    -   *Linux/WSL*: Uses apt (`sudo apt install hmmer`).
2.  Check/Install **Java 17+** (required for Java PPRFinder compatibility).
3.  Check **Python 3.10+** and create a virtual environment.
4.  Install all dependencies.

### 2. Launch App

**Mac**:
- Double-click `start.command`.

**Linux / Windows (WSL)**:
- Open terminal and run:
  ```bash
  bash start.command
  ```

The application will open automatically in your default web browser.

## Quick Start (Docker)

Build image:

```bash
docker build -t ppprpp:latest .
```

Run container:

```bash
docker run --rm -it \
  -p 7860:7860 \
  -v "$(pwd)/data/input:/app/data/input" \
  -v "$(pwd)/data/output:/app/data/output" \
  ppprpp:latest
```

Then open: `http://127.0.0.1:7860`

## Troubleshooting

- **"Python 3.10 or newer is required"**: Install Python 3.10+ from [python.org](https://www.python.org/downloads/) and retry.
- **"hmmsearch not found"**: The installer attempts to install HMMER. If it fails, please install it manually:
    -   **Mac**: `brew install hmmer`
    -   **Linux/WSL**: `sudo apt install hmmer`
- **"Java 17+ not found"**: Install Java 17 and retry.
    -   **Mac**: `brew install openjdk@17`
    -   **Linux/WSL**: `sudo apt install openjdk-17-jre-headless`
    -   Optional: set `PPRFINDER_JAVA` to the Java 17 binary path (e.g. `/opt/homebrew/opt/openjdk@17/bin/java`).
- **Browser doesn't open**: Look for a URL in the console output (e.g., `http://127.0.0.1:7860`) and open it manually.
- **PS_Scan returns no output or fails**: The PS_Scan option uses the ExPASy ScanProsite REST API (https://prosite.expasy.org), which requires an active internet connection. If the API is unavailable (network error, server downtime), an error message will appear in the System Log at the bottom of the interface. In this case, switch to **PPRFinder** (the default, recommended tool) which runs locally without network access.

## Basic Usage

### 1. Inputs (Tab 1: Inputs)
This tab is for preparing your data.

1.  **Protein Input**: 
    - Enter the Protein Sequence (FASTA format) or upload a file. This sequence is used for Motif Discovery.
2.  **Input Mode**: Select how you want to provide the Target RNA.
    - **Manual Sequence**: Paste an RNA sequence directly.
    - **NCBI Accession**: Enter an NCBI ID (e.g., `AP017386`) to automatically fetch the FASTA and Annotation (GFF3).
3.  **Custom Annotations** (Optional):
    - You can manually add gene information (Name, Type, Start, End, Strand) to be displayed in the analysis.
    - Click **"Save (GFF3)"** to save your custom genes.
4.  **Preview & Features**:
    - **Preview Genes**: Visualize the features (Genes, CDS, etc.) in the specified range.
    - **Range**: Enter Start/End to zoom in. Leave empty or 0 to view the whole sequence (in Manual Mode).
    - **Feature Types**: Select which feature types to display (e.g., `gene`, `CDS`, `other`).

### 2. Analysis (Tab 2: Settings & Run)
Configure the analysis parameters.

1.  **Discovery Tool**: Choose the motif discovery algorithm.
    - **PPRFinder** (default, recommended): Runs locally using HMMER. No internet connection required.
    - **PS_Scan**: Uses the ExPASy ScanProsite API. Requires an active internet connection; if the server is unavailable, an error is shown in the System Log.
2.  **Scoring Method**: Select **P-type** or **PLS-type** scoring. Choose based on the class of your PPR protein.
3.  **Parameters**:
    - **Gap Size**: Set the exact gap size to use (e.g., 3). The analysis evaluates candidate windows matching this gap size, together with Gap 0 as a baseline. A value of 0 means no gap (contiguous motif matching only).
    - **Use Motif Confidence Score (+Score)**: Enable to weight scores by per-motif confidence values from the discovery tool.
4.  **Visual Settings**:
    - **Plot Height**: Adjust the height of the result genome browser.
    - **Highlight Top**: Highlights the top-scoring binding site with a star marker.
    - **Show Strand**: Filter the browser display. Select **Forward**, **Reverse**, or **Both**. When **Both** is selected, forward and reverse strand scores are shown in separate subplots. The active filter is shown in the plot title (e.g., "Strand: + only"). Note: if the target RNA is single-stranded or only one strand is scored, switching to the other strand will show an empty plot.
5.  **Run**:
    - Enter a Job Name (optional).
    - Click **"Run Analysis"**.

### 3. Results (Tab 3: Results)
View the analysis output.

- **Status Log**: Shows progress and key statistics.
- **PPR Motif Plot**: Visualizes the discovered PPR motifs and their amino acids.
- **Genome Browser**:
    - Interactive plot showing binding scores across the RNA sequence.
    - Displays gene annotations (from NCBI or Custom).
    - Click on any point to inspect detailed binding info.
- **Detailed Inspection**:
    - Shows the specific alignment between the PPR motifs and the RNA sequence at the clicked position.
- **Downloads**:
    - Download the results as CSV, BedGraph, BED, or Reference FASTA/GFF.

### 4. File Management
- Use the **"📂 Input Folder"** and **"📂 Output Folder"** buttons at the top right to quickly open the data directories on your computer.

## Benchmark Evaluation (CLI)

For reproducible evaluation (off-target, S/N, gap sweep), run:

```bash
python scripts/benchmark_prediction.py \
  --protein-fasta /path/to/protein.fasta \
  --rna-fasta /path/to/rna.fasta \
  --known-sites /path/to/known_sites.bed \
  --tool PPRFinder \
  --method yan \
  --max-gap 3 \
  --output-dir data/benchmark
```

Notes:
- `PPRFinder` uses the Java implementation (`PPRfinder.jar`) for reference-compatible motif discovery.
- `--known-sites-coordinate` supports `auto`, `bed`, `zero-based-half-open`, and `one-based-inclusive`.
- `--gap-values 0,1,2,3` can be used instead of `--max-gap`.

Main outputs:
- `data/benchmark/benchmark_summary.csv`: site-level metrics (validated score, rank percentile, S/N, off-target counts).
- `data/benchmark/gap_sweep.csv`: gap-level summary (runtime, detected rate, median S/N, etc.).
- `data/benchmark/benchmark_run_meta.json`: run configuration for reproducibility.

## Citation

If you use PPPRPP in your research, please cite:

> Igarashi, K. (2026). PPPRPP: A web application for P-type PPR protein target site prediction with gap-tolerant scoring. *Plant Biotechnology*, PB-S-1341.

The source code is archived on Zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18371990.svg)](https://doi.org/10.5281/zenodo.18371990)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
