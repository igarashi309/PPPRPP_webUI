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

def plot_genes_only(gff_path, start, end, title_id="Sequence", feature_type="gene"):
    """
    Plots only the gene annotations within a range using Plotly.
    """
    from . import ncbi
    import plotly.graph_objects as go
    import os
    
    fig = go.Figure()
    
    # Base line
    fig.add_trace(go.Scatter(
        x=[start, end], y=[0, 0], 
        mode='lines', line=dict(color='gray', width=1), 
        showlegend=False, hoverinfo='skip'
    ))
    
    if not gff_path or not os.path.exists(gff_path):
        fig.update_layout(title=f"No GFF file found for {title_id}")
        return fig
        
    # Parse GFF
    # Assuming ncbi.parse_gff is available and works
    try:
        features = ncbi.parse_gff(gff_path, start, end)
    except Exception as e:
        fig.update_layout(title=f"Error parsing GFF: {e}")
        return fig
        
    # Filter features
    # Support single string or list
    types_to_show = feature_type if isinstance(feature_type, list) else [feature_type]
    
    # Standard types usually have dedicated colors
    STANDARD_TYPES = ["gene", "CDS", "tRNA", "rRNA", "exon"]
    
    # Constants for plotting
    gene_y_fwd = 0.5
    gene_y_rev = -0.5
    thickness = 0.4
    half_thick = thickness / 2
    
    # Arrow head size
    view_width = end - start
    if view_width <= 0: view_width = 100
    head_len = view_width * 0.015 # 1.5% of view
    
    annotations = []
    
    for feat in features:
        f_type = feat['type']
        
        # Filtering Logic
        show_it = False
        if f_type in types_to_show:
            show_it = True
        elif "other" in types_to_show:
            if f_type not in STANDARD_TYPES:
                show_it = True
        
        if not show_it:
            continue
            
        f_start = max(start, feat['start'])
        f_end = min(end, feat['end'])
        
        # Skip if out of view (parse_gff might return overlaps)
        if f_end < start or f_start > end:
            continue
            
        f_strand = feat['strand']
        
        # Name
        attrs = feat.get('attributes', {})
        f_name = attrs.get('Name', attrs.get('gene', attrs.get('ID', f_type)))
        
        # Color
        if f_type == 'CDS': color = 'orange'
        elif f_type == 'tRNA': color = 'green'
        elif f_type == 'rRNA': color = 'purple'
        elif f_type == 'gene': color = 'blue'
        elif f_type == 'exon': color = 'cornflowerblue'
        else: color = 'gray' # Default for 'other'
        
        y_pos = gene_y_fwd if f_strand == '+' else gene_y_rev
        
        # Draw Arrow Shape using Scatter Fill
        # Points: Start,Top -> HeadBase,Top -> Tip,Mid -> HeadBase,Bot -> Start,Bot
        
        # Adjust lengths for aesthetic arrow
        draw_start = f_start
        draw_end = f_end
        
        # Ensure minimum visibility
        if draw_end - draw_start < view_width * 0.001:
             draw_end = draw_start + view_width * 0.001
        
        h_len = min(head_len, draw_end - draw_start)
        
        y_top = y_pos + half_thick
        y_bot = y_pos - half_thick
        y_mid = y_pos
        
        if f_strand == '+':
            x_tip = draw_end
            x_base = draw_end - h_len
            x_tail = draw_start
            
            x_pts = [x_tail, x_base, x_tip, x_base, x_tail, x_tail]
            y_pts = [y_top, y_top, y_mid, y_bot, y_bot, y_top]
        else:
            x_tip = draw_start
            x_base = draw_start + h_len
            x_tail = draw_end
            
            x_pts = [x_tail, x_base, x_tip, x_base, x_tail, x_tail]
            y_pts = [y_top, y_top, y_mid, y_bot, y_bot, y_top]
            
        fig.add_trace(go.Scatter(
            x=x_pts, y=y_pts,
            fill="toself",
            mode='lines',
            line=dict(width=1, color='black'),
            fillcolor=color,
            name=f"{f_name} ({f_strand})",
            hoverinfo='text',
            hovertext=f"{f_name}<br>{f_type}<br>{f_start}-{f_end} ({f_strand})"
        ))
        
        # Label (Centered)
        label_x = (draw_start + draw_end) / 2
        annotations.append(dict(
            x=label_x, y=y_pos,
            text=f_name,
            showarrow=False,
            font=dict(color='white' if f_type in ['gene', 'CDS'] else 'black', size=10),
            bgcolor=None # Transparent background inside colored arrow
        ))
        
    fig.update_layout(
        title=f"Sequence Features: {title_id} ({start}-{end})",
        yaxis=dict(
            range=[-1.5, 1.5], 
            showgrid=False, 
            zeroline=False, 
            showticklabels=False,
            title="Strand (+ / -)"
        ),
        xaxis=dict(title="Position", range=[start, end]),
        height=300,
        annotations=annotations,
        showlegend=False
    )
    
    return fig
