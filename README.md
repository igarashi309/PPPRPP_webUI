# P-type PPR Prediction Program (PPPRPP)

A standalone web application for identifying PPR motifs in protein sequences and analyzing their binding sites on RNA targets.
For `PPRFinder`, this app uses the Java implementation (`PPRfinder.jar`) so motif discovery matches the reference behavior.

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
2.  Check/Install **Java 17+** (required for Java PPRFinder compatibility).
3.  Create a Python environment.
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

- **"hmmsearch not found"**: The installer attempts to install HMMER. If it fails, please install it manually:
    -   **Mac**: `brew install hmmer`
    -   **Linux/WSL**: `sudo apt install hmmer`
- **"Java 17+ not found"**: Install Java 17 and retry.
    -   **Mac**: `brew install openjdk@17`
    -   **Linux/WSL**: `sudo apt install openjdk-17-jre-headless`
    -   Optional: set `PPRFINDER_JAVA` to the Java 17 binary path (e.g. `/opt/homebrew/opt/openjdk@17/bin/java`).
- **Browser doesn't open**: Look for a URL in the console output (e.g., `http://127.0.0.1:7860`) and open it manually.

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

1.  **Tool Selection**: Choose the analysis algorithm (e.g., `Kobayashi 2000`, `Yan 2019`).
2.  **Parameters**:
    - **Gap Size**: Set the maximum gap size to consider (e.g., 3).
    - **Scoring**: Enable/Disable P-S Score usage.
3.  **Visual Settings**:
    - **Plot Height**: Adjust the height of the result genome browser.
    - **Highlight Scores**: Option to highlight the top-scoring binding site.
4.  **Run**:
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

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
