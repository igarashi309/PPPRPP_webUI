# How to Visualize Results in IGV

The PPR Motif Analyzer generates files compatible with the **Integrative Genomics Viewer (IGV)**. This allows you to explore binding sites, scores, and gene annotations interactively.

## 1. Install IGV
If you haven't already, download and install the IGV Desktop application from the official website:
[https://igv.org/doc/desktop/](https://igv.org/doc/desktop/)

## 2. Prepare Your Files
From the "NCBI Analysis" tab in the PPR Motif Analyzer, download the following files:
1.  **FASTA**: The reference genome sequence (e.g., `AP012527.fasta`).
2.  **BedGraph**: The calculated binding scores (e.g., `AP012527_kobayashi_gap0_ps.bedgraph`).
3.  **BED**: The binding sites visualized as reads (e.g., `AP012527_kobayashi_gap0_ps.bed`).
4.  **GFF3**: The gene annotations (e.g., `AP012527.gff3`).

## 3. Load Files into IGV

### Step 1: Load the Genome (FASTA)
1.  Open IGV.
2.  Go to **Genomes** > **Load Genome from File...**
3.  Select the downloaded **FASTA file** (`.fasta`).
    *   *Note: Do not use "File > Load from File" for the FASTA. It must be loaded as a genome.*

### Step 2: Load Data Tracks (BedGraph, BED, & GFF3)
1.  Go to **File** > **Load from File...**
2.  Select the **BedGraph**, **BED**, and **GFF3** files. You can select all at once.

## 4. Explore the Data
- **Navigation**: Use the search bar to zoom into specific regions or genes.
- **Tracks**:
    - **Reference**: Shows the DNA sequence (zoom in to see bases).
    - **Annotations (GFF3)**: Shows gene locations, CDS, tRNAs, etc.
    - **PPR Binding Score (BedGraph)**: Shows the binding score as a bar chart. Higher bars indicate stronger predicted binding.
    - **PPR Binding Sites (BED)**: Shows binding sites as blocks.
        - **Blue**: Forward strand binding.
        - **Red**: Reverse strand binding.
        - **Hover**: Hover over a block to see the PPR name and score.
- **Data Range**: If the BedGraph bars look too small or cut off:
    - Right-click the BedGraph track.
    - Select **Autoscale** to fit the data to the view.
    - Or select **Set Data Range** to manually set the min/max (e.g., 0 to 100).

## Tips
- **Color**: The BedGraph track is blue by default. You can change it by right-clicking the track > **Change Track Color**.
- **Overlay**: You can load multiple BedGraph files (e.g., different gap settings) to compare them.
