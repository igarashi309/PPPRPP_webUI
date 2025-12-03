"""
Visualization functions for PPR analysis.
"""
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import logomaker
import plotly.graph_objects as go
import pandas as pd
import numpy as np
from .config import KOBAYASHI_BINDING_SCORES, YAN_BINDING_SCORES

def plot_motif_sequence_logo(motif_list, binding_scores_dict, title="Sequence Logo (Bits)", fig_width=8, fig_height=2):
    """
    Plots a sequence logo for a list of motifs.
    """
    data = []
    for motif in motif_list:
        if motif in binding_scores_dict:
            p = np.array(binding_scores_dict[motif], dtype=np.float64).tolist()
        else:
            p = np.array([0.25, 0.25, 0.25, 0.25], dtype=np.float64).tolist()
        data.append(p)
    df_logo = pd.DataFrame(data, columns=["A", "C", "G", "U"])

    def calc_entropy(arr):
        arr = np.array(arr, dtype=np.float64)
        mask = arr > 0
        return -np.sum(arr[mask] * np.log2(arr[mask]))

    df_info = df_logo.copy()
    for i in df_logo.index:
        H = calc_entropy(df_logo.loc[i])
        I = 2 - H
        df_info.loc[i] = df_logo.loc[i] * I

    logo = logomaker.Logo(df_info, shade_below=0.5, fade_below=0.5, font_name='Arial Rounded MT Bold')
    logo.ax.set_title(title)
    logo.ax.set_xlabel("Motif Position")
    logo.ax.set_ylabel("Information (bits)")
    logo.ax.set_xticks(np.arange(len(motif_list)))
    logo.ax.set_xticklabels(motif_list, rotation=90)
    logo.fig.set_size_inches(fig_width, fig_height)
    logo.fig.tight_layout()
    plt.show()

def plot_binding_scores(df, title, highlight_positions, height=400):
    """
    Plots binding scores across the sequence.
    """
    df['PPR_Binding_Score'] = df['PPR_Binding_Score'].astype(float)
    max_value = df['PPR_Binding_Score'].max() * 1.1
    
    unique_ppr_names = df['PPR_Name'].unique()
    # Updated to use matplotlib.colormaps
    cmap = plt.get_cmap('coolwarm')
    colors = [cmap(i) for i in np.linspace(0, 1, len(unique_ppr_names))]

    plt.figure(figsize=(15, 5))
    
    for i, ppr_name in enumerate(unique_ppr_names):
        ppr_df = df[df['PPR_Name'] == ppr_name]
        plt.scatter(ppr_df['Position'], ppr_df['PPR_Binding_Score'], color=colors[i], label=ppr_name, alpha=0.5, s=10)

    top1_index = df['PPR_Binding_Score'].idxmax()
    plt.scatter(df['Position'][top1_index], df['PPR_Binding_Score'][top1_index],
                color='yellow', edgecolor='black', s=50, label='Top 1 Score', marker='*')

    for hp in highlight_positions:
        plt.axvline(x=hp, color='purple', linestyle='--', label='Highlight')

    plt.title(title)
    plt.xlabel("Position")
    plt.ylabel("PPR Binding Score")
    plt.ylim(0, max_value)
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
    plt.tight_layout()
    plt.show()

def plot_top_and_worst_scores(df, gap_info, motif_info, binding_source, height=400):
    """
    Plots detailed scores for the top and worst binding sites using Plotly.
    """
    df['PPR_Binding_Score'] = df['PPR_Binding_Score'].astype(float)
    
    if motif_info == "Motif score included":
        # Handle string representation of lists in CSV
        if isinstance(df['PPR_Motif_Scores'].iloc[0], str):
             df['PPR_Motif_Scores_List'] = df['PPR_Motif_Scores'].apply(lambda x: list(map(float, x.split(','))))
        else:
             df['PPR_Motif_Scores_List'] = df['PPR_Motif_Scores']
             
        all_scores = df['PPR_Motif_Scores_List'].explode()
        max_value = all_scores.max() * 1.1
    else:
        max_value = 1

    max_binding_score_row = df.loc[df['PPR_Binding_Score'].idxmax()]
    min_binding_score_row = df.loc[df['PPR_Binding_Score'].idxmin()]

    if motif_info == "Motif score included":
        scores_max = max_binding_score_row['PPR_Motif_Scores_List']
        scores_min = min_binding_score_row['PPR_Motif_Scores_List']
    else:
        # Handle string representation
        if isinstance(max_binding_score_row['PPR_Binding_Scores_List'], str):
            scores_max = [float(s) for s in max_binding_score_row['PPR_Binding_Scores_List'].split(',')]
            scores_min = [float(s) for s in min_binding_score_row['PPR_Binding_Scores_List'].split(',')]
        else:
            scores_max = max_binding_score_row['PPR_Binding_Scores_List']
            scores_min = min_binding_score_row['PPR_Binding_Scores_List']

    motif_labels_gap_max = max_binding_score_row['Motif_Labels_GAP'].split(',')
    segment_max = max_binding_score_row['Segment']
    motif_labels_gap_min = min_binding_score_row['Motif_Labels_GAP'].split(',')
    segment_min = min_binding_score_row['Segment']

    def prepare_labels(motif_labels):
        motif_count = 0
        labels_with_indices = []
        numbers = []
        for label in motif_labels:
            if label == "--":
                labels_with_indices.append("")
                numbers.append("")
            else:
                motif_count += 1
                labels_with_indices.append(label)
                numbers.append(f"{motif_count}")
        return numbers, labels_with_indices

    mn_max, ml_max = prepare_labels(motif_labels_gap_max)
    mn_min, ml_min = prepare_labels(motif_labels_gap_min)

    # Plotly Top
    fig_max = go.Figure()
    fig_max.add_trace(go.Bar(
        x=[f"{mn}<br>{seg}<br>{ml}" for mn, seg, ml in zip(mn_max, segment_max, ml_max)],
        y=scores_max,
        text=[f"{s:.2f}" for s in scores_max],
        texttemplate='%{text}',
        textposition='auto',
        marker_color='indianred'
    ))
    fig_max.update_layout(
        title=f'【Top Scores】{binding_source}, {gap_info}, {motif_info}',
        yaxis=dict(range=[0, max_value]),
        height=height
    )
    fig_max.show()

    # Plotly Worst
    fig_min = go.Figure()
    fig_min.add_trace(go.Bar(
        x=[f"{mn}<br>{seg}<br>{ml}" for mn, seg, ml in zip(mn_min, segment_min, ml_min)],
        y=scores_min,
        text=[f"{s:.2f}" for s in scores_min],
        texttemplate='%{text}',
        textposition='auto',
        marker_color='dodgerblue'
    ))
    fig_min.update_layout(
        title=f'【Worst Scores】{binding_source}, {gap_info}, {motif_info}',
        yaxis=dict(range=[0, max_value]),
        height=height
    )
    fig_min.show()

def plot_top_rna_sequence_logomaker(top_segment, title="Top Score RNA Sequence", fig_width=10, fig_height=2):
    """
    Plots RNA sequence logo.
    """
    seq_rna = top_segment.upper().replace("T", "U").strip()
    letters = list(seq_rna)
    
    rows = []
    for letter in letters:
        row = {"A": 0, "C": 0, "G": 0, "U": 0}
        if letter in row:
            row[letter] = 2
        rows.append(row)
    df_rna = pd.DataFrame(rows, columns=["A", "C", "G", "U"])

    logo = logomaker.Logo(df_rna, shade_below=0, fade_below=0, font_name='Arial Rounded MT Bold', color_scheme='classic')
    logo.ax.set_title(title)
    logo.ax.set_xlabel("Position")
    logo.ax.set_xticks(np.arange(len(letters)))
    logo.ax.set_xticklabels(np.arange(1, len(letters)+1), rotation=90)
    logo.fig.set_size_inches(fig_width, fig_height)
    logo.fig.tight_layout()
    plt.show()

def plot_combined_top_logos(top_motif_labels, binding_dict_for_logo, top_segment, binding_source, gap_info, motif_info, fig_width=10, fig_height=4):
    """
    Combines bits logo and RNA sequence logo.
    """
    # Top Panel: Bits Logo
    data = []
    for motif in top_motif_labels:
        if motif in binding_dict_for_logo:
            p = np.array(binding_dict_for_logo[motif], dtype=np.float64)
        else:
            p = np.array([0.25, 0.25, 0.25, 0.25], dtype=np.float64)
        data.append(p)
    df_logo = pd.DataFrame(data, columns=["A", "C", "G", "U"])

    def calc_entropy(arr):
        arr = np.array(arr, dtype=np.float64)
        mask = arr > 0
        return -np.sum(arr[mask] * np.log2(arr[mask]))

    df_info = df_logo.copy()
    for i in df_logo.index:
        H = calc_entropy(df_logo.loc[i])
        I = 2 - H
        df_info.loc[i] = df_logo.loc[i] * I

    # Bottom Panel: RNA Sequence
    seq_rna = top_segment.upper().replace("T", "U").strip()
    letters = list(seq_rna)
    rna_rows = []
    for letter in letters:
        row = {"A": 0, "C": 0, "G": 0, "U": 0}
        if letter in row:
            row[letter] = 2
        rna_rows.append(row)
    df_rna = pd.DataFrame(rna_rows, columns=["A", "C", "G", "U"])

    fig, (ax_top, ax_bottom) = plt.subplots(nrows=2, ncols=1, figsize=(fig_width, fig_height), sharex=False)

    logo_top = logomaker.Logo(df_info, ax=ax_top, shade_below=0.5, fade_below=0.5, font_name='Arial Rounded MT Bold')
    ax_top.set_title(f"Bits Sequence Logo (Top Score) | {binding_source}, {gap_info}, {motif_info}")
    ax_top.set_ylabel("Information (bits)")
    ax_top.set_xticks(np.arange(len(top_motif_labels)))
    ax_top.set_xticklabels(top_motif_labels)

    logo_bottom = logomaker.Logo(df_rna, ax=ax_bottom, shade_below=0.5, fade_below=0.5, font_name='Arial Rounded MT Bold', color_scheme='classic')
    ax_bottom.set_title("Top Score RNA Sequence")
    ax_bottom.set_xlabel("Position")
    ax_bottom.set_ylabel("One-hot value")
    ax_bottom.set_xticks(np.arange(len(letters)))
    ax_bottom.set_xticklabels(np.arange(1, len(letters)+1))

    fig.subplots_adjust(hspace=0.5)
    fig.suptitle(f"Combined Top Score Logos", fontsize=12)
    plt.tight_layout(rect=[0, 0, 1, 1])
    plt.show()
