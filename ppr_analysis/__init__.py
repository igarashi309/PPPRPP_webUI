"""
PPR Analysis Package
"""
from .config import KOBAYASHI_BINDING_SCORES, YAN_BINDING_SCORES
from .core import analyze_ppr_data, process_fasta_sequences

try:
    from .io import read_scores_from_file, load_ppr_data_from_pickle, save_scores_to_csv
except ImportError:
    print("Warning: Could not import IO module (requires pandas).")

try:
    from .visualization import (
        plot_motif_sequence_logo, 
        plot_binding_scores, 
        plot_top_and_worst_scores, 
        plot_top_rna_sequence_logomaker, 
        plot_combined_top_logos
    )
except ImportError:
    print("Warning: Could not import Visualization module (requires matplotlib, logomaker, plotly).")

try:
    from .utils import infer_additional_info, get_top_n_per_ppr
except ImportError:
    print("Warning: Could not import Utils module (requires pandas).")
