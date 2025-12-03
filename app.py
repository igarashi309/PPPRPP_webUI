import gradio as gr
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import os
import tempfile
import pickle
import numpy as np
from pathlib import Path

# Import package modules
from ppr_analysis import core, io as ppr_io, visualization, config, motif_generation, ncbi

# Check Gradio Version
if int(gr.__version__.split('.')[0]) < 4:
    raise RuntimeError(f"Gradio version {gr.__version__} is too old. Please install gradio>=4.0.0")
print(f"Gradio Version: {gr.__version__}")

# --- Logic Functions ---

def run_unified_analysis(protein_text, protein_file, rna_text, rna_file, tool_choice, analysis_type, job_name):
    # 1. Process Protein Input (Discovery)
    if protein_file is not None:
        with open(protein_file.name, 'r') as f:
            protein_input = f.read()
    else:
        protein_input = protein_text

    if not protein_input:
        return "Error: Please provide a protein sequence.", None, None, None, None, None

    if not protein_input.strip().startswith(">"):
        protein_input = f">UserSequence\n{protein_input}"

    # 2. Process RNA Input (Analysis Target)
    if rna_file is not None:
        with open(rna_file.name, 'r') as f:
            rna_input = f.read()
    else:
        rna_input = rna_text

    if not rna_input:
        return "Error: Please provide an RNA sequence.", None, None, None, None, None
        
    sequence = ""
    seq_id = job_name
    if rna_input.strip().startswith(">"):
        parts = rna_input.strip().split("\n", 1)
        sequence = parts[1].replace("\n", "").strip()
        # If job name is default, try to use fasta header
        if job_name == "my_analysis":
            seq_id = parts[0][1:].strip()
    else:
        sequence = rna_input.strip()

    # 3. Run Discovery
    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".fasta") as tmp_fasta:
        tmp_fasta.write(protein_input)
        tmp_fasta_path = tmp_fasta.name

    output_dir = tempfile.mkdtemp()
    pkl_path = None
    status_msg = ""

    try:
        if tool_choice == "PPRFinder":
            pkl_path = motif_generation.run_pprfinder(tmp_fasta_path, output_dir)
            # Motifs file is now saved as {seq_id}_motifs.txt in output_dir
            seq_id = os.path.basename(pkl_path).replace('.pkl', '')
            motifs_file = Path(output_dir) / f"{seq_id}_motifs.txt"
        else: # PS_Scan
            # PS_Scan log parsing is not fully implemented in this version to produce .pkl
            # We will try to run it but might fail to proceed to analysis if .pkl is not generated.
            # For unified flow, we strongly recommend PPRFinder.
            log_path = motif_generation.run_ps_scan(tmp_fasta_path, output_dir)
            return f"PS_Scan executed. Log saved to {os.path.basename(log_path)}. (Unified analysis requires PPRFinder for now)", None, None, None, None, None, None, None
            
        if not pkl_path:
             return "Error: Motif discovery failed.", None, None, None, None, None, None, None

    except Exception as e:
        return f"Error during discovery: {e}", None, None, None, None, None, None, None
    finally:
        if os.path.exists(tmp_fasta_path):
            os.remove(tmp_fasta_path)

    # 4. Run Analysis
    status_msg += f"\nRunning Analysis ({analysis_type})..."
    
    # Load the discovered motifs
    try:
        motif_data = ppr_io.load_ppr_data_from_pickle(pkl_path)
        ppr_name = os.path.basename(pkl_path)
    except Exception as e:
        return f"Error loading motifs: {e}", None, None, None, None, None, None, None

    # Run Analysis (Both Strands)
    min_front = 7
    min_back = 7
    max_gap = 0
    
    # Prepare RNA sequence
    # rna_input is already defined from rna_text or rna_file earlier.
    if not rna_input:
        return "Error: No RNA sequence provided.", None, None, None, None, None, None, None

    # Parse RNA FASTA
    try:
        rna_record = next(SeqIO.parse(io.StringIO(rna_input), "fasta"))
        rna_seq = str(rna_record.seq).upper().replace("U", "T") # Convert to DNA for analysis if needed, or keep RNA. Core usually expects DNA/RNA string.
    except Exception as e:
        # Fallback if raw sequence or malformed FASTA
        rna_seq = rna_input.strip().upper().replace("U", "T")

    all_scores = []
    predicted_scores = [] # Default empty, or could allow manual input. For unified, we assume default/empty.
    
    try:
        # Forward
        scores_fwd = core.analyze_ppr_data(motif_data, rna_seq, min_front, min_back, max_gap, predicted_scores, ppr_name)
        if analysis_type in scores_fwd:
            for s in scores_fwd[analysis_type]:
                row = list(s)
                row.append("Forward")
                row.append(row[0]) # Genome_Position for forward is same as relative
                all_scores.append(row)
        
        # Reverse
        rev_seq = core.get_reverse_complement(rna_seq)
        scores_rev = core.analyze_ppr_data(motif_data, rev_seq, min_front, min_back, max_gap, predicted_scores, ppr_name)
        if analysis_type in scores_rev:
            for s in scores_rev[analysis_type]:
                row = list(s)
                row.append("Reverse")
                genome_pos = len(rna_seq) - 1 - row[0] # Adjust position for reverse strand
                row.append(genome_pos)
                all_scores.append(row)
                
    except Exception as e:
        import traceback
        traceback.print_exc()
        return f"Error during analysis: {e}", None, None, None, None, None, None, None

    if not all_scores:
        return "Error: No binding sites found.", None, None, None, None, None, None, None

    # Create DataFrame
    columns = [
        "Position_Relative", "PPR_Binding_Score", "PPR_Name", "Motif_Labels",
        "Motif_Labels_GAP", "Gap_Size", "Segment", "PPR_Motif_Scores",
        "Total_Score", "PPR_Binding_Scores_List", "Predicted PPR Motif Scores List",
        "Strand", "Genome_Position"
    ]
    df = pd.DataFrame(all_scores, columns=columns)

    # CSV
    csv_path = f"{seq_id}_{analysis_type}.csv"
    df.to_csv(csv_path, index=False)

    # 5. Visualization (Genome Browser)
    fig = plot_genome_browser(df, seq_id, analysis_type)
    
    # 6. Visualize Motifs
    motif_df, motif_fig = pd.DataFrame(), go.Figure()
    if motifs_file and os.path.exists(motifs_file):
        motif_df, motif_fig = visualize_motifs(motifs_file)
    else:
        status_msg += "\nWarning: Motif details file not found for visualization."

    return "Analysis Complete", df, csv_path, fig, go.Figure(), df, motif_df, motif_fig
    
def plot_genome_browser(df, seq_id="Sequence", analysis_type="Analysis"):
    """Helper to plot the genome browser."""
    fig = go.Figure()
    
    unique_ppr_names = df['PPR_Name'].unique()
    colors = px.colors.qualitative.Plotly
    
    for i, name in enumerate(unique_ppr_names):
        color = colors[i % len(colors)]
        
        # Forward
        df_fwd = df[(df['PPR_Name'] == name) & (df['Strand'] == "Forward")]
        if not df_fwd.empty:
            fig.add_trace(go.Scatter(
                x=df_fwd['Genome_Position'],
                y=df_fwd['PPR_Binding_Score'],
                mode='markers',
                name=f"{name} (+)",
                marker=dict(color=color, size=6, symbol='circle'),
                customdata=df_fwd.index, # Store index for click event
                hovertemplate="Pos: %{x}<br>Score: %{y:.2f}<br>Strand: +<extra></extra>"
            ))
            
        # Reverse
        df_rev = df[(df['PPR_Name'] == name) & (df['Strand'] == "Reverse")]
        if not df_rev.empty:
            fig.add_trace(go.Scatter(
                x=df_rev['Genome_Position'],
                y=-df_rev['PPR_Binding_Score'],
                mode='markers',
                name=f"{name} (-)",
                marker=dict(color=color, size=6, symbol='x'),
                customdata=df_rev.index,
                hovertemplate="Pos: %{x}<br>Score: %{y:.2f}<br>Strand: -<extra></extra>"
            ))

    fig.update_layout(
        title=f"Genome Browser - {seq_id}",
        xaxis_title="Genomic Position",
        yaxis_title="Binding Score (+ / -)",
        height=600,
        hovermode="closest"
    )
    return fig

def visualize_motifs(motifs_file):
    """Parses motifs file and returns DataFrame and Plotly figure."""
    try:
        # Read _motifs.txt (created by ppr_finder.py)
        # Columns: header, start, end, score, motif_seq, second, fifth, last, name
        # No header in file
        col_names = ['Header', 'Start', 'End', 'Score', 'Sequence', '2nd', '5th', 'Last', 'Type']
        df = pd.read_csv(motifs_file, sep='\t', names=col_names, header=None)
        
        if df.empty:
             return pd.DataFrame(), go.Figure()

        # Add Code column
        df['Code'] = df['5th'] + df['Last']
        
        # Plot
        fig = go.Figure()
        
        # Draw protein backbone (line)
        max_pos = df['End'].max() + 50
        fig.add_trace(go.Scatter(
            x=[0, max_pos], y=[0, 0],
            mode='lines',
            line=dict(color='gray', width=2),
            showlegend=False,
            hoverinfo='skip'
        ))
        
        # Color map
        color_map = {
            'P': '#1f77b4', 'P1': '#aec7e8', 'P2': '#ff7f0e',
            'L': '#2ca02c', 'L1': '#98df8a', 'L2': '#d62728',
            'S': '#9467bd', 'S1': '#c5b0d5', 'S2': '#8c564b',
            'SS': '#e377c2'
        }
        
        # Draw motifs as rectangles
        for _, row in df.iterrows():
            color = color_map.get(row['Type'], '#bcbd22')
            
            # Rectangle
            fig.add_shape(
                type="rect",
                x0=row['Start'], x1=row['End'],
                y0=-0.2, y1=0.2,
                line=dict(color="black", width=1),
                fillcolor=color,
                opacity=0.8
            )
            
            # Text (Code)
            fig.add_annotation(
                x=(row['Start'] + row['End']) / 2,
                y=0,
                text=row['Code'],
                showarrow=False,
                font=dict(color="white", size=10)
            )
            
            # Text (Type) - above
            fig.add_annotation(
                x=(row['Start'] + row['End']) / 2,
                y=0.3,
                text=row['Type'],
                showarrow=False,
                font=dict(color="black", size=8)
            )

        fig.update_layout(
            title="Predicted PPR Motifs",
            xaxis_title="Amino Acid Position",
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, range=[-0.5, 0.5]),
            height=250,
            margin=dict(l=20, r=20, t=40, b=20)
        )
        
        return df[['Type', 'Start', 'End', 'Score', '5th', 'Last', 'Code', 'Sequence']], fig
        
    except Exception as e:
        print(f"Error visualizing motifs: {e}")
        return pd.DataFrame(), go.Figure()

def generate_detail_figure(row):
    """Helper to generate the detail bar chart for a single row."""
    try:
        def parse_list_col(val):
            if isinstance(val, list): return val
            return [float(x) for x in str(val).split(',')]

        scores = parse_list_col(row['PPR_Motif_Scores'])
        motif_labels = row['Motif_Labels_GAP'].split(',') if isinstance(row['Motif_Labels_GAP'], str) else row['Motif_Labels_GAP']
        segment = row['Segment']
        
        # Prepare data for plotting
        x_labels = []
        colors = []
        text_labels = []
        
        # Define color map for common PPR motifs
        color_map = {
            'P': '#1f77b4',  # Blue
            'P1': '#aec7e8', # Light Blue
            'P2': '#ff7f0e', # Orange
            'L': '#2ca02c',  # Green
            'L1': '#98df8a', # Light Green
            'L2': '#d62728', # Red
            'S': '#9467bd',  # Purple
            'S1': '#c5b0d5', # Light Purple
            'S2': '#8c564b', # Brown
            'SS': '#e377c2', # Pink
            '--': '#7f7f7f'  # Gray (Gap)
        }
        
        for i, label in enumerate(motif_labels):
            # X-axis: Nucleotide + Index (1-based relative to site)
            nuc = segment[i] if i < len(segment) else "-"
            x_labels.append(f"{i+1}<br><b>{nuc}</b>")
            
            # Color
            colors.append(color_map.get(label, '#bcbd22')) # Default to Olive if unknown
            
            # Text on bar
            text_labels.append(label)

        fig = go.Figure(go.Bar(
            x=x_labels, 
            y=scores, 
            text=text_labels,
            textposition='auto', 
            marker_color=colors
        ))
        
        fig.update_layout(
            title=f"Detail: {row['PPR_Name']} at {row['Genome_Position']} ({row['Strand']})",
            xaxis_title="Position / Nucleotide",
            yaxis_title="Binding Score",
            height=400,
            showlegend=False
        )
        return fig
    except Exception as e:
        return go.Figure(layout=dict(title=f"Error generating plot: {e}"))

def update_detail_plot(evt: gr.SelectData, df):
    if df is None or df.empty:
        return go.Figure()
        
    try:
        x_val = evt.value['x']
        y_val = evt.value['y']
        
        mask = (df['Genome_Position'] == x_val) & \
               ((np.isclose(df['PPR_Binding_Score'], y_val)) | (np.isclose(df['PPR_Binding_Score'], -y_val)))
        
        row = df[mask]
        
        if row.empty:
            return go.Figure(layout=dict(title="No data found for clicked point"))
            
        return generate_detail_figure(row.iloc[0])
        
    except Exception as e:
        return go.Figure(layout=dict(title=f"Error: {e}"))

def show_manual_detail(position, df):
    """Show detail for a manually entered position."""
    if df is None or df.empty:
        return go.Figure(layout=dict(title="No analysis data available"))
    
    if position is None:
        return go.Figure(layout=dict(title="Please enter a position"))

    try:
        # Find rows at this position
        # We might have multiple (different PPRs). We'll pick the one with max score.
        rows = df[df['Genome_Position'] == int(position)]
        
        if rows.empty:
            return go.Figure(layout=dict(title=f"No binding sites found at position {position}"))
            
        # Sort by score descending and pick top
        best_row = rows.sort_values('PPR_Binding_Score', ascending=False).iloc[0]
        
        return generate_detail_figure(best_row)
    except Exception as e:
        return go.Figure(layout=dict(title=f"Error: {e}"))

# --- UI Definition ---

with gr.Blocks(title="PPR Motif Analyzer (Unified)") as demo:
    gr.Markdown("# PPR Motif Analyzer (Unified Workflow)")
    gr.Markdown("Input your Protein sequence (to discover motifs) and RNA sequence (target) to analyze binding sites.")
    
    with gr.Row():
        with gr.Column():
            gr.Markdown("### 1. Protein Input (Motif Discovery)")
            prot_text = gr.Textbox(label="Protein FASTA", lines=5, placeholder=">ProteinID\\nSEQUENCE...")
            prot_file = gr.File(label="Or Upload Protein File")
            tool_select = gr.Dropdown(choices=["PPRFinder", "PS_Scan"], value="PPRFinder", label="Discovery Tool")
            
        with gr.Column():
            gr.Markdown("### 2. RNA Input (Target Analysis)")
            rna_text = gr.Textbox(label="RNA FASTA", lines=5, placeholder=">RNA_ID\\nSEQUENCE...")
            rna_file = gr.File(label="Or Upload RNA File")
            
    with gr.Row():
        with gr.Column():
            gr.Markdown("### 3. Settings & Run")
            analysis_opts = [
                'kobayashi_gap0_ps', 'kobayashi_gap_ps', 'kobayashi_gap0_no_ps', 'kobayashi_gap_no_ps',
                'yan_gap0_ps', 'yan_gap_ps', 'yan_gap0_no_ps', 'yan_gap_no_ps'
            ]
            analysis_select = gr.Dropdown(choices=analysis_opts, value='kobayashi_gap0_ps', label="Analysis Type")
            job_name_input = gr.Textbox(value="my_analysis", label="Job Name")
            run_btn = gr.Button("Run Unified Analysis", variant="primary")
            
        with gr.Column():
            status_output = gr.Textbox(label="Status")
            
    # Outputs
    with gr.Row():
        with gr.Column():
            gr.Markdown("### Motif Discovery Results")
            motif_plot = gr.Plot(label="Motif Architecture")
            motif_table = gr.DataFrame(label="Motif Details")
            
    with gr.Row():
        plot_browser = gr.Plot(label="Genome Browser (Click point to see details)")
    
    with gr.Row():
        with gr.Column(scale=1):
            gr.Markdown("### Inspect Specific Site")
            manual_pos = gr.Number(label="Genome Position", precision=0)
            manual_btn = gr.Button("Show Detail")
        with gr.Column(scale=3):
            plot_detail = gr.Plot(label="Selected Site Detail")
        
    with gr.Accordion("Data Table & Download", open=False):
        df_output = gr.DataFrame(label="Results Summary")
        csv_output = gr.File(label="Download CSV")
        
    # State
    df_state = gr.State()
    
    # Events
    run_btn.click(
        run_unified_analysis,
        inputs=[prot_text, prot_file, rna_text, rna_file, tool_select, analysis_select, job_name_input],
        outputs=[status_output, df_output, csv_output, plot_browser, plot_detail, df_state, motif_table, motif_plot]
    )
    
    print(f"DEBUG: plot_browser type: {type(plot_browser)}")
    print(f"DEBUG: plot_browser attributes: {dir(plot_browser)}")
    
    try:
        plot_browser.select(update_detail_plot, inputs=[df_state], outputs=[plot_detail])
    except AttributeError:
        print("WARNING: plot_browser.select() failed. Interactive features disabled.")
    except Exception as e:
        print(f"WARNING: Failed to attach select event: {e}")
        
    # Manual Event
    manual_btn.click(
        show_manual_detail,
        inputs=[manual_pos, df_state],
        outputs=[plot_detail]
    )

if __name__ == "__main__":
    demo.launch(share=True)
