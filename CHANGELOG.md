# Changelog

## [v0.0.4] - 2025-12-05

### Added
- **Range Inspector**: New feature in NCBI Analysis tab to visualize specific genomic ranges.
    - Displays binding scores and gene annotations.
    - Symmetric Y-axis for easier comparison of + and - strands.
    - Gene annotations displayed as arrows to indicate direction.
    - Adjustable plot height (400px - 1200px).
    - Feature filtering (gene, CDS, tRNA, rRNA, exon).
- **IGV Support**:
    - Generated files (BedGraph, BED, GFF3) are optimized for IGV.
    - `IGV_INSTRUCTIONS.md` added.
    - BED file generation for "read-like" visualization of binding sites.
- **Temp File Cleanup**: Automatically cleans up generated analysis files (*.csv, *.fasta, etc.) on startup.

### Changed
- **NCBI Analysis UI**: Simplified to remove heavy DataFrame and Genome Browser components, preventing browser crashes.
- **Y-Axis Scaling**: Unified Y-axis scaling across all Genome Browser plots (symmetric limits).
- **GFF3 Handling**: Flattened GFF3 files (removed Parent attributes from exons/CDS) for better visualization in IGV.

### Fixed
- **Performance**: Optimized core analysis loop using NumPy for significant speedup.
- **Stability**: Fixed Gradio timeouts by enabling queue and truncating UI data tables.
- **Bugs**:
    - Fixed `NameError` in `run_unified_analysis`.
    - Fixed `AttributeError` in `run_ncbi_analysis`.
    - Fixed FASTA parsing position shift issue.
    - Fixed "Inspect Specific Site" error when handling empty lists.
    - Removed problematic Range Slider to prevent `TypeError`.
- **Security**: Disabled public URL (share link) by default.
- **Visualization**: Added `Analysis_Type` to "Inspect Specific Site" plot title.
