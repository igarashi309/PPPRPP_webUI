# Changelog

## [v1.1.1] - 2026-07-16

### Added
- **Vector (SVG) Figure Export**: The "Download plot" (camera) button on all Plotly figures now exports **SVG** (vector) by default instead of PNG, giving publication-quality output that scales without pixelation and is editable in Illustrator/Inkscape. Implemented centrally by injecting a Plotly `config` into Gradio's plot serialization — no new dependencies and no changes to individual figure code. (Client-side export supports SVG/PNG/JPEG/WebP; PDF is not available from the button — convert the SVG if a PDF is required.)

### Fixed
- **Detail Inspector Axis Alignment**: In the Detail Inspector plot, the left (Binding Score) and right (Motif Confidence) y-axes auto-ranged independently, placing their zero baselines at different vertical positions. This made bars appear as if Binding Score could be negative. Both axes now use `rangemode='tozero'`, anchoring the zero of Motif Confidence to the zero of Binding Score. Binding Score values themselves are unchanged (always non-negative); this is a display-only fix.

## [v1.1.0] - 2026-03-04

### Added
- **Benchmark Evaluation Script**: Added `scripts/benchmark_prediction.py` for gap-size sweep benchmarking. Computes S/N ratio (validated-site score / Q95 of off-target scores), rank percentile, detected rate, and runtime across user-specified gap sizes.
- **Revision Figure Generator**: Added `scripts/generate_revision_figures.py` to produce Supplementary Figures S2/S3 and Supplementary Table S1 from benchmark CSV outputs.

### Changed
- **Dependency Upgrade**: Upgraded `gradio` from `4.44.1` to `6.2.0`. Updated `demo.launch()` arguments and component APIs to align with Gradio 6.x specifications.
- **Gap Size Parameter Label**: Unified user-facing label from `max_gap` / `Allow Gaps` to **Gap Size** throughout the UI and documentation for consistency with the manuscript.

### Fixed
- **PPRFinder Scoring Bug**: Corrected a motif-score calculation error in the PPRFinder integration. Previously, the combined score was incorrectly computed; the fix ensures accurate per-window scores used for target-site ranking.

## [v0.1.1] - 2026-01-20

### Fixed
- **App Stability**: Pinned `gradio` version to `4.44.1` to prevent compatibility issues with newer updates (v6.x).
- **Launch Error**: Corrected `demo.launch()` arguments to align with Gradio 4.x specifications, fixing startup crashes on some environments.

## [v0.1.0] - 2026-01-19

### Added
- **Interactive Gap Selection**: In "Inspect Specific Site", users can now see all alignment variations (different gap sizes and positions) for a specific site. The dropdown allows switching between them instantly.
- **RNA Sequence Support**: All inputs (Sequence Text, FASTA files) are now automatically normalized to Uppercase RNA (U-notation), allowing flexible input of DNA/RNA mixed sequences.
- **RNA Display**: Detail Inspector and plots now consistently display sequences using RNA notation (U).

### Changed
- **UX Improvement**: "Alignment Variation" dropdown is now always visible, making the feature more discoverable.
- **Robustness**: Fixed data type issues preventing "Inspect Specific Site" from working in Lite Mode context.
- **Stability**: Fixed a startup crash caused by incorrect State initialization order.

## [v0.0.8] - 2026-01-09

### Changed
- **Result Filenames**: filenames now include detailed analysis parameters (e.g., `JobName_ProtID_RNAID_P-type_Gap0_ScoreOn_results.csv`) for better organization and identification in "Past Results".
- **Motif Architecture Plot**: 
    - Title now displays analysis settings (Method, Gap, Score).
    - Confidence Score legend moved to top-right and horizontal for better visibility.

## [v0.0.6] - 2025-12-05

### Fixed
- **NCBI Analysis**: Fixed a critical bug where NCBI Analysis failed due to incorrect function signature when calling `run_unified_analysis` internally.

## [v0.0.5] - 2025-12-05

### Changed
- **Analysis Selection UI**: Replaced the single dropdown with intuitive Radio buttons and Checkboxes for "Scoring Method", "Allow Gaps", and "Use Motif Confidence Score". This allows for 8 clear combinations of analysis parameters.

### Fixed
- **Score Calculation**: Fixed a bug where `Total_Score` incorrectly included the Motif Confidence Score (PS Score) even when "Use Motif Confidence Score" was disabled.

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
