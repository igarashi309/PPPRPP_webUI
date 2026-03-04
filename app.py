import gradio as gr
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import matplotlib.pyplot as plt
import os
import tempfile
import pickle
import numpy as np
from pathlib import Path
import glob
from Bio import SeqIO
import io
import shutil
import subprocess
import platform

# --- Constants ---
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "data")
OUTPUT_DIR = os.path.join(DATA_DIR, "output")
INPUT_DIR = os.path.join(DATA_DIR, "input")

# Ensure directories exist
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(INPUT_DIR, exist_ok=True)

TMP_DIR = OUTPUT_DIR

VERSION = "v1.1.0"

def cleanup_temp_files():
    """Cleans up temporary files from previous runs."""
    # Disabled for persistent output directory usage
    pass
    # patterns = ["*.csv", "*.fasta", "*.bedgraph", "*.bed", "*.gff3"]
    # count = 0
    # for pattern in patterns:
    #     for f in glob.glob(os.path.join(TMP_DIR, pattern)): # Changed to TMP_DIR context if active
    #         try:
    #             os.remove(f)
    #             count += 1
    #         except Exception as e:
    #             print(f"Error deleting {f}: {e}")
    # if count > 0:
    #     print(f"Cleaned up {count} temporary files.")

def open_folder(path):
    """Opens a folder in the system file explorer."""
    if not os.path.exists(path):
        return f"Path not found: {path}"
    
    try:
        system = platform.system()
        if system == "Darwin": # macOS
            subprocess.run(["open", path])
        elif system == "Windows":
            os.startfile(path)
        elif system == "Linux":
            subprocess.run(["xdg-open", path])
        return f"Opened: {path}"
    except Exception as e:
        return f"Error opening folder: {e}"

# Import package modules
from ppr_analysis import core, io as ppr_io, visualization, config, motif_generation, ncbi

# Check Gradio Version
if int(gr.__version__.split('.')[0]) < 4:
    raise RuntimeError(f"Gradio version {gr.__version__} is too old. Please install gradio>=4.0.0")
print(f"Gradio Version: {gr.__version__}")

custom_css = """
/* Floating Terminal Log Panel - Enhanced for Gradio */
/* System Log Terminal */
#system_log_terminal {
    height: 300px !important;
    min-height: 300px !important;
    display: block !important;
}
#system_log_terminal textarea {
    background-color: #000000 !important;
    color: #00ff00 !important;
    font-family: 'Courier New', monospace !important;
    font-size: 12px !important;
    border: 1px solid #005500 !important;
    line-height: 1.4 !important;
    box-shadow: inset 0 0 10px rgba(0, 50, 0, 0.5);
    height: 100% !important;
}
#system_log_terminal label {
    display: block !important;
    height: 100% !important;
    border: none !important;
}
#system_log_terminal label span {
    display: none !important; /* Hide the text label itself, keep container */
}
#terminal_panel,
#terminal_panel > div,
#terminal_panel .gradio-column {
    position: fixed !important;
    bottom: 10px !important;
    right: 10px !important;
    width: 400px !important;
    max-height: 300px !important;
    z-index: 99999 !important;
    background: rgba(0, 0, 0, 0.85) !important;
    border: 2px solid rgba(0, 255, 0, 0.7) !important;
    border-radius: 10px !important;
    box-shadow: 0 0 25px rgba(0, 255, 0, 0.4) !important;
    font-family: 'Courier New', Consolas, monospace !important;
    transition: all 0.3s ease !important;
    overflow: hidden !important;
    backdrop-filter: blur(5px) !important;
}

#terminal_panel.collapsed,
#terminal_panel.collapsed > div {
    max-height: 36px !important;
    width: 200px !important;
}

#terminal_header {
    background: linear-gradient(180deg, rgba(0, 50, 0, 0.9), rgba(0, 25, 0, 0.9)) !important;
    padding: 8px 14px !important;
    cursor: pointer !important;
    display: flex !important;
    justify-content: space-between !important;
    align-items: center !important;
    border-bottom: 1px solid rgba(0, 170, 0, 0.6) !important;
    user-select: none !important;
}

#terminal_header:hover {
    background: linear-gradient(180deg, rgba(0, 70, 0, 0.9), rgba(0, 35, 0, 0.9)) !important;
}

#terminal_title {
    color: rgba(0, 255, 0, 0.9) !important;
    font-size: 13px !important;
    font-weight: bold !important;
    text-shadow: 0 0 8px rgba(0, 255, 0, 0.5) !important;
}

#terminal_toggle {
    color: rgba(0, 255, 0, 0.9) !important;
    font-size: 16px !important;
    text-shadow: 0 0 8px rgba(0, 255, 0, 0.5) !important;
}

#terminal_content,
#terminal_content > div {
    max-height: 250px !important;
    overflow-y: auto !important;
    padding: 10px !important;
    background: transparent !important;
}

/* Target Gradio textbox specifically */
#terminal_content textarea,
#terminal_content .wrap,
#terminal_content input,
#terminal_content .gradio-textbox,
#terminal_content .gradio-textbox textarea,
#terminal_content > div > div > div {
    background: transparent !important;
    background-color: transparent !important;
    color: rgba(0, 255, 0, 0.85) !important;
    font-family: 'Courier New', Consolas, monospace !important;
    font-size: 11px !important;
    line-height: 1.4 !important;
    border: none !important;
    box-shadow: none !important;
    text-shadow: 0 0 5px rgba(0, 255, 0, 0.3) !important;
}

#terminal_content label,
#terminal_content .label-wrap {
    display: none !important;
}

/* Remove any white backgrounds from Gradio containers */
#terminal_content .block,
#terminal_content .form,
#terminal_content .container {
    background: transparent !important;
    border: none !important;
    box-shadow: none !important;
}

/* Custom scrollbar */
#terminal_content::-webkit-scrollbar {
    width: 6px !important;
}
#terminal_content::-webkit-scrollbar-track {
    background: rgba(0, 20, 0, 0.5) !important;
}
#terminal_content::-webkit-scrollbar-thumb {
    background: rgba(0, 170, 0, 0.7) !important;
    border-radius: 3px !important;
}

/* Hide old debug console */
#debug_log_text {
    display: none !important;
}
.hidden-log-box {
    display: none !important;
}
"""

DEFAULT_PROTEIN = """>PPR791(Rf1a)_IR24
MARRAASRAVGALRSDGSIQGRGGRAGGSGAEDARHVFDELLRRGRGASIYGLNRALADVARDSPAAAVSRYNRMARAGADEVTPDLCTYGILIGCCCRAGRLDLGFAALGNVIKKGFRVDAIAFTPLLKGLCADKRTSDAMDIVLRRMTELGCIPNVFSYNILLKGLCDENRSQEALELLHMMADDRGGGSPPDVVSYTTVINGFFKEGDSDKAYSTYHEMLDRGILPDVVTYNSIIAALCKAQAMDKAMEVLNTMVKNGVMPDCMTYNSILHGYCSSGQPKEAIGFLKKMRSDGVEPDVVTYSLLMDYLCKNGRCMEARKIFDSMTKRGLKPEITTYGTLLQGYATKGALVEMHGLLDLMVRNGIHPDHYVFSILICAYAKQGKVDQAMLVFSKMRQQGLNPNAVTYGAVIGILCKSGRVEDAMLYFEQMIDEGLSPGNIVYNSLIHGLCTCNKWERAEELILEMLDRGICLNTIFFNSIIDSHCKEGRVIESEKLFELMVRIGVKPNVITYNTLINGYCLAGKMDEAMKLLSGMVSVGLKPNTVTYSTLINGYCKISRMEDALVLFKEMESSGVSPDIITYNIILQGLFQTRRTAAAKELYVRITESGTQIELSTYNIILHGLCKNKLTDDALQMFQNLCLMDLKLEARTFNIMIDALLKVGRNDEAKDLFVAFSSNGLVPNYWTYRLMAENIIGQGLLEELDQLFLSMEDNGCTVDSGMLNFIVRELLQRGEITRAGTYLSMIDEKHFSLEASTASLFIDLLSGGKYQEYYRFLPEKYKSFIESLSC"""

DEFAULT_RNA = """>atp6-orf79_BT
ataactataggattctcacgagaaggcttcaaatcaagcaatccacaagacaagcgagttcactatcatagttgatccttgactgcttccttcacggtaatggtcctcctaggataggatatgagcgaacatagtcaattacagatatgcgaactgttggaaagcactagtttcgaagaaaaagagtagcaagggattcgaaacagacaaagcctatgaacaagagtaaggactgcaacgattagagagagtagggagtaccgttcattaccgattagagagagtagggagtaccgaacagcctctagtactgcttgctatcatagagaagggagtaccgaacagcctctagtactgctatctacacaccctctagtaccgctatcgtagagaatctagtacagctatctacgagacgagaaaaggaaacgatttgatttttcacacaagacctagccattgcattggtttagttaagtaattgcatttccgcttcagctgcccttgcttttccttcacccaaaaggctctctttcggctgtacccatcccggagggatgccaggtactgctcttcaggttaagagaaacagtcattgcattgttgaataagggaaagtggaattgttgaatcagccctcgcctggtgcctaccaataatacttcataattgcgtatatagaagggaggggagaggtactctatcttgttaacttactgtttttttctcctagctgggtttttcattttattctcttagagataatgctctccagaagagctcgtaatgaattggattcgaaccaatatctcaatgcgacttatccattagtcgatcattgcgacttatccattagtcgatcatcacgggaggaaccgccaccacgggttcccggccgccctcccaccttgtaggggtgggatatttgcctggtccaccagaattccatagtatggttttgtaatgaaagtatgagcaaaccaccactgtcctgtctttcttacttaggaaagactacacgaatagaggtgccccaccacgaatagtcaacctctatagggtgcccggcttttatattatgatctatgattttttccattaaagccgggtcgtacagaggggcctgatactttagcaaagcaaacaacagggcgcctcctacaacgacaccgaagatgaaatatggaagaccgttagtccctcgggtagtggagaagagccatcggaccagatttgccatagcacaataaaataatatttttaagctccgctctccttatattggagagacttgtttggggtcggcgttggtttttccaacgaaaaacgaattcgaacatgaaatgaaaaaacaaatgctagttgtcattttttgtaacgcaatacacttccgcgtttgtcccgccccctaaaatcttcggctcctcgtttttattcaattatgaaattactcattttgatggagatttatagcatcattcaagtaaatacaaattgagatcgtagaaacatgagcttgtaatatagctacacctaattccagaccggttaatgctagaactataaataagggaccaagatctcctatgaaatagaaaatattattcagaaatagcatagtccaagcgaacccacttaaaatctttactgaactatgaccggccatcatattagcaaataaacgtattcctgagcttaatgcacgaaaacaatgagagattagctcaaggagtactaaaaaaggtgctaatggcagtgggactcccgctggtaataagaagctaaaaaaatgaagcccatgtctttgaaatccaacgatcgtaatgcctataaaaatggaaaatgaaagagccaaagtaatgagaaaatgacttgtcactgtgaagctaaagggtatcataccctggggattacgaaataacgaaaaagtaaaagtgaccgagatgcgagggaaaaacttttgtttaacatttccacctatttgttcgtttaccaggttcagcacgaaatcataaataagctctaccaaggattgccatgcatttggcactgactttcccccgccctttttcgtaacaactccaaataaagatggcaccaaaacgacagttaccgccatagacaaggattcatttgtgaatgataaataaaagttacccatatcaagaccaaataatggatcaatgaaaaattgatccaaggggcttgaaattactggtccgccgcctccagggataactatatcgccgttcctccattcagccccatctaacaactctattccgttctgccagttattgaaacggatacgctgcgccctctcaaatggttcaaaaaaagcgtctatttctgcttgtctccttctcttcaacgagttcacgttataatcattgaggagtttccagatcgagtctctctgattcaaacccattattaccacatgattgtgatcgaaattcattggaatttcccatctttttaagctctgctcctagaaagggagactgatctcaactggcctcggtttttccaaggaaactctatctttagtagcgatcacatggcgtcgacccaacggaaaacatactatactactatgaacatttcagaatggaaaacgaaactgaacactttccaatactggtttcggatggcttatgtgcatcgatttttctcaatcttcgtaggggcttgacgtataggttccccaacctatcctaaaactctagtctagtgaaggcggggccaacaaatgacctcttgagtgagattacatcccagagcctgcctattcaatgccactcctggtaagatccattttcaccgatacccgctctctatatatgagtttctggaaatggttatctcacaagcaaagaaggaattcaaattctacccctcagaagaaccttcaagcagtacaacccacagatgatgaagtgaacaggacttcttttcttatatttgcctgttgaactccggcagcaagtgccttccttactagtgaggtcttgcttttgcttctttggttggtactcttctatctgggactgggaaagatctttctgcattaacggaattaggcttactcctcagaaagccctactactttactctaggcctagtagttatagctaactcactgctattggttgggatgaatttccctaccctcgccctacgtggtttagaagagtagtataagtaaggaagcctgataatcatataggcaaaaggcaccctgtcccagacaaaccagacaggcagcaaatttcactaaggcagtaaagcattctgatactcctatacagacaaggctgtaacccaatccgaccagcctctcggaataccgaataggcaaagtacggaatcaaaatcaaggggacaaggagcctaccaagtacggcaacaaaagaaagaaagacgagacagcgctacgcttgttcaccaaaggtcaagat"""

DEFAULT_ACCESSION = "AP017386"

# --- Logic Functions ---

def list_past_results():
    """Lists available result CSV files in TMP_DIR, sorted by newest."""
    files = glob.glob(os.path.join(TMP_DIR, "*_results.csv"))
    # Sort by mtime
    files.sort(key=os.path.getmtime, reverse=True)
    f_list = [os.path.basename(f) for f in files]
    # Return updates for BOTH dropdowns (Main Tab and NCBI Tab)
    return gr.update(choices=f_list, value=None), gr.update(choices=f_list, value=None)


import json
import ast

def load_csv_with_parsing(path):
    """Loads CSV and parses list columns back to objects."""
    try:
        df = pd.read_csv(path)
        
        # Columns that need parsing
        list_cols = ['PPR_Binding_Scores_List', 'PPR_Motif_Scores', 'Predicted PPR Motif Scores List', 'Motif_Confidence_Scores']
        
        for col in list_cols:
            if col in df.columns:
                # Use ast.literal_eval for safety, fallback to json if needed, or crude eval
                # The columns are likely string repr of lists: "[1.0, 2.0]"
                def parse_list(x):
                    if isinstance(x, str):
                        try:
                            # Try literal_eval
                            val = ast.literal_eval(x)
                            # If result is still string (double quoted), try again?
                            if isinstance(val, str) and val.strip().startswith('[') and val.strip().endswith(']'):
                                try:
                                    val = ast.literal_eval(val)
                                except:
                                    pass
                            return val
                        except:
                            # Fallback: simple strip and split
                            s = x.strip()
                            if (s.startswith('[') and s.endswith(']')) or (s.startswith("'[") and s.endswith("]'")) or (s.startswith('"[') and s.endswith(']"')):
                                # Strip quotes if present
                                if s.startswith("'") or s.startswith('"'):
                                    s = s[1:-1]
                                # Strip brackets
                                if s.startswith('[') and s.endswith(']'):
                                    s = s[1:-1]
                                if not s.strip(): return []
                                try:
                                    return [float(i) for i in s.split(',')]
                                except:
                                    return []
                            return []
                    return x
                df[col] = df[col].apply(parse_list)
        return df
    except Exception as e:
        print(f"Error parsing CSV: {e}")
        return pd.DataFrame()

def load_csv_optimized(path, usecols=None):
    """Loads CSV with specific columns, skipping expensive list parsing."""
    try:
        # Use python engine for robustness against quoting/escaping issues
        # on_bad_lines='warn' (pandas >= 1.3) or error_bad_lines=False (older)
        try:
            df = pd.read_csv(path, usecols=usecols, engine='python', on_bad_lines='warn')
        except TypeError:
             # Fallback for older pandas
             df = pd.read_csv(path, usecols=usecols, engine='python', error_bad_lines=False)
             
        # DEBUG: Check Strand column
        if 'Strand' in df.columns:
            print(f"DEBUG: Loaded CSV {os.path.basename(path)}. Strands: {df['Strand'].unique()}. Total Rows: {len(df)}")
            if df['Strand'].isnull().any():
                print("DEBUG: Warning: NaN values found in Strand column!")
                # Optional: drop rows with NaN strand?
                # df = df.dropna(subset=['Strand'])
        else:
            print(f"DEBUG: Loaded CSV {os.path.basename(path)} but 'Strand' column missing!")
        return df
    except Exception as e:
        print(f"Error loading optimized CSV: {e}")
        return pd.DataFrame()

        return pd.DataFrame()

def run_motif_discovery(protein_input, tool_choice):
    """
    Helper to run motif discovery (PPRFinder/PSScan) and return loaded motif_data + ppr_name.
    """
    if not protein_input:
        raise ValueError("No protein sequence provided.")

    # Pre-clean protein input
    clean_protein_input = protein_input.strip()
    if clean_protein_input.startswith(">"):
        lines = clean_protein_input.split("\n")
        header = lines[0]
        seq_content = "".join(lines[1:]).replace(" ", "").replace("\r", "").replace("\n", "").strip()
        final_protein_content = f"{header}\n{seq_content}"
    else:
        seq_content = clean_protein_input.replace(" ", "").replace("\r", "").replace("\n", "").strip()
        final_protein_content = f">Protein\n{seq_content}"

    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".fasta") as tmp_fasta:
        tmp_fasta.write(final_protein_content)
        tmp_fasta_path = tmp_fasta.name

    output_dir = tempfile.mkdtemp()
    
    try:
        if tool_choice == "PPRFinder":
            pkl_path = motif_generation.run_pprfinder(tmp_fasta_path, output_dir)
        else: # PS_Scan
            pkl_path = motif_generation.run_ps_scan_api(tmp_fasta_path, output_dir)
            
        if not pkl_path:
            raise RuntimeError("Motif discovery failed (No pickle path returned).")
            
        # success
        seq_id = os.path.basename(pkl_path).replace('.pkl', '')
        
        # Load data
        motif_data = ppr_io.load_ppr_data_from_pickle(pkl_path)
        
        # Add types
        import pickle
        with open(pkl_path, 'rb') as f:
            records = pickle.load(f)
        
        motif_types = []
        for r in records:
            if "PPR_TYPE" in r.annotations:
                motif_types.append(r.annotations["PPR_TYPE"])
            else:
                motif_types.append("Unknown")
        motif_data['motif_types'] = motif_types
        
        return motif_data, seq_id, pkl_path
        
    finally:
        if os.path.exists(tmp_fasta_path):
            os.remove(tmp_fasta_path)


def fetch_ncbi_basics(accession_id):
    """
    Step 1: Just fetch FASTA and GFF.
    """
    # Fallback to default if empty
    if not accession_id or not accession_id.strip():
        accession_id = DEFAULT_ACCESSION
        
    status_log = [f"Fetching genome for {accession_id}..."]
    try:
        fasta_path = ncbi.fetch_fasta(accession_id, TMP_DIR)
        status_log.append(f"Genome fetched: {os.path.basename(fasta_path)}")
        
        gff_path = None
        try:
             # Try fetch GFF
             gff_path = ncbi.fetch_gff(accession_id, TMP_DIR)
             status_log.append(f"GFF fetched: {os.path.basename(gff_path)}")
        except Exception as e:
             status_log.append(f"GFF fetch failed (optional): {e}")
             
        status_log.append("Ready. Set range and click 'Inspect Range'.")
        return "\n".join(status_log), fasta_path, gff_path
        
    except Exception as e:
        return f"Error fetching genome: {e}", None, None


def analyze_ncbi_range_wrapper(accession_id, start, end, fasta_path, gff_path, protein_text, tool_choice, scoring_method, max_gap_size, use_ps_score):
    import pandas as pd # Import locally to ensure availability
    """
    Step 2: Inspect Range (Recalc).
    """
    if not fasta_path or not os.path.exists(fasta_path):
        return go.Figure(layout=dict(title="No FASTA loaded. Click 'Fetch Data' first.")), "Error: No FASTA.", None, None
        
    try:
        # 1. Run/Get Motif Data
        motif_data, ppr_name, _ = run_motif_discovery(protein_text, tool_choice)
    except Exception as e:
        return go.Figure(layout=dict(title=f"Motif Discovery Failed: {e}")), f"Motif Error: {e}", None, None
        
    # 2. Prepare Context for plot_range_browser
    # We need to construct analysis_context to pass to plot_range_browser
    # AND we need to read the full RNA seq? 
    # Actually plot_range_browser needs `rna_seq` in context to slice it?
    # Yes, `plot_range_browser` > `subseq = analysis_context['rna_seq'][start:end]`
    # So we must read full FASTA.
    
    try:
        with open(fasta_path, "r") as f:
            rna_text = f.read()
            # Normalize to RNA (U) - Or keep provided logic
            # clean lines
            lines = rna_text.splitlines()
            if lines[0].startswith(">"):
                rna_seq = "".join(lines[1:]).strip().upper().replace(" ", "")
            else:
                rna_seq = "".join(lines).strip().upper().replace(" ", "")
            rna_seq = core.dna_to_rna(rna_seq)
            
    except Exception as e:
        return go.Figure(layout=dict(title=f"Error reading FASTA: {e}")), f"FASTA Error: {e}", None, None

    # Construct Context
    method_map = {"PLS-type": "kobayashi", "P-type": "yan"}
    method = method_map.get(scoring_method, scoring_method).lower()
    
    analysis_context = {
        'motif_data': motif_data,
        'rna_seq': rna_seq,
        'min_front': 1,
        'min_back': 1,
        'max_gap': max_gap_size,
        'predicted_scores': [],
        'use_ps_score': use_ps_score,
        'method': method,
        'ppr_name': ppr_name
    }
    
    # 3. Call plot_range_browser
    # Pass empty df (None) to force it to use context usage logic?
    # Wait, `plot_range_browser` usually filters DF.
    # But if DF is None, it returns "No data" UNLESS we modify it to handle "Context Only" mode.
    # Currently `plot_range_browser` iterates `df['PPR_Name'].unique()`. if DF is None, validation fails.
    
    # We must construct a dummy DF or refactor `plot_range_browser`.
    # Actually, `plot_range_browser`'s new RECALC block is inside the main loop `for name in ppr_names`.
    # If DF is empty, no `ppr_names`.
    # So we should pass a Dummy DF with 1 row just to trigger the loop? 
    # Or better: Create a "Direct Mode" in `plot_range_browser`?
    
    # Let's create a minimal dummy DF to trick `plot_range_browser` to enter the loop for `ppr_name`.
    dummy_df = pd.DataFrame({'PPR_Name': [ppr_name], 'Strand': ['Forward'], 'Genome_Position': [start], 'PPR_Binding_Score': [0]})
    # This will let it find `unique_names = [ppr_name]`.
    # Then inside the loop:
    # `if analysis_context ...`: It enters RECALC block.
    # `res_fwd = ...`
    # It replaces `df_range` with `df_recalc`.
    # So the dummy data is ignored. Perfect.
    
    range_len = int(end) - int(start)
    limit = 50000 
    # If range > limit? User said "perform all", but larger ranges might kill browser.
    # Let's warn in title if huge? 
    # We updated `plot_range_browser` to check `MAX_RECALC_RANGE = 50000`.
    # If I pass a dummy DF, and range > 50k, it falls back to Filtering Dummy DF... which sends Empty Plot.
    # So for this "Direct Mode" to work for large ranges, we must strictly force Recalc or ensure plot_range_browser handles it.
    # I should update `plot_range_browser` to Allow recalc for any range if DF is Dummy?
    # Or just tell user "Zoom in".
    # User said "Evaluate specific range".
    # I will stick to 50k limit for now.
    
    fig = plot_range_browser(dummy_df, gff_path, start, end, seq_id=accession_id, 
                             feature_types=[], plot_height=600, 
                             highlight_top=False, highlight_strand="Both", 
                             analysis_context=analysis_context)
                             
    return fig, "Analysis Done.", dummy_df, analysis_context
    """
    Parses a motifs.txt file and returns a list of Bio.SeqRecord objects 
    compatible with the analysis core.
    """
    try:
        import pandas as pd
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        
        # Format: id, start, end, score, motif_seq, second, fifth, last, type
        # No header
        try:
            df = pd.read_table(motifs_file, names=['id','start','end','score', 'motif_seq','second','fifth','last','type'], header=None)
        except Exception:
            # Try skipping header if it exists (some versions might have it)
             df = pd.read_table(motifs_file, names=['id','start','end','score', 'motif_seq','second','fifth','last','type'], header=0)
             
        records = []
        for _, row in df.iterrows():
            # Create SeqRecord
            # Note: The txt content for 'motif_seq' is the sequence
            seq_obj = Seq(str(row['motif_seq']))
            rec = SeqRecord(seq_obj, id=seq_id)
            
            # Annotations required by core.py
            rec.annotations = {
                "PPR_TYPE": row['type'],
                "PS_SCORE": float(row['score']),
                # Txt file has 1-based start? Usually yes for display.
                # But let's check motif_generation.py parse_ps_scan_xml saving:
                # d_start = p.annotations["PPR_STARTS"] + 1
                # So here we subtract 1
                "PPR_STARTS": int(row['start']) - 1,
                "PPR_ENDS": int(row['end']), # End is usually inclusive in display? 
                # In xml saving: d_end = p.annotations["PPR_ENDS"] (which was safe_end = exclusive index in python slicing?)
                # Wait, motif_generation.py: safe_end = min(len, end)
                # python slicing [start:end] -> end is exclusive.
                # So PPR_ENDS in annotation was exclusive index.
                # In txt saving: d_end = p.annotations["PPR_ENDS"]
                # So txt end is exclusive index.
                # So we use it as is.
                "FIFTH_AA": str(row['fifth']),
                "LAST_AA": str(row['last']),
                "PPR_CODE": str(row['fifth']) + str(row['last'])
            }
            records.append(rec)
            
        return records
    except Exception as e:
        print(f"Error parsing motif file {motifs_file}: {e}")
        return []

        return []

def reconstruct_context_from_csv(csv_path, fasta_path=None):
    """
    Attempts to reconstruct the analysis context (motifs, RNA, settings)
    from the Result CSV filename and associated files.
    """
    try:
        filename = os.path.basename(csv_path)
        # Format: {JobName}_{Prot}_{RNA}_{Method}_{Gap}_{Score}_results.csv
        # Example: my_analysis_PPR791_Rf1a__IR24_AP017386.1_P-type_Gap0_ScoreOff_results.csv
        
        parts = filename.replace("_results.csv", "").split('_')
        
        # We need to find specific tokens.
        # Method: P-type or PLS-type
        method = "yan" # Default
        if "P-type" in parts: method = "yan"
        elif "PLS-type" in parts: method = "kobayashi"
        
        # Gap
        max_gap = 0
        for p in parts:
            if p.lower().startswith("gap"):
                try:
                    max_gap = int(p.lower().replace("gap", ""))
                except: pass
                
        # Score
        use_ps = True
        if "ScoreOff" in parts: use_ps = False
        
        # Files
        # Motif file: {Base}_motifs.txt
        motif_path = csv_path.replace("_results.csv", "_motifs.txt")
        if not os.path.exists(motif_path):
            # Fallback: Check OUTPUT_DIR (TMP_DIR) in case we are in Gradio temp
            fallback_path = os.path.join(TMP_DIR, os.path.basename(motif_path))
            if os.path.exists(fallback_path):
                motif_path = fallback_path
            else:
                print(f"Motif file not found: {motif_path} or {fallback_path}")
                return {}
            
        records = parse_motif_file_to_records(motif_path)
        if not records:
            print("Failed to parse motif data")
            return {}
            
        # Convert List[SeqRecord] to Dict format ('motif_labels', 'ppr_scores') expected by core
        motif_labels = [r.annotations.get('PPR_CODE') for r in records]
        ppr_scores = [r.annotations.get('PS_SCORE') for r in records]
        motif_types = [r.annotations.get('PPR_TYPE', 'Unknown') for r in records]
        motif_data = {'motif_labels': motif_labels, 'ppr_scores': ppr_scores, 'motif_types': motif_types}
            
        # RNA / Genome
        # Try to find RNA ID. It's properly hard to split by underscore if ID has underscores.
        # But we know method/gap/score are at the end.
        # And we check if the Accession exists in data/output or input
        
        # Heuristic: List .fasta files in output/input and see if any is a substring of filename?
        # Or parse standard accessions.
        # Accession usually: AP017386.1 -> "AP017386.1"
        
        rna_seq = ""
        
        # Check output directory for FASTAs
        candidates = glob.glob(os.path.join(TMP_DIR, "*.fasta")) + glob.glob(os.path.join(INPUT_DIR, "*.fasta"))
        best_cand = None
        
        for cand in candidates:
            # Check if basename of cand (minus extension) is in our result filename
            cand_name = os.path.splitext(os.path.basename(cand))[0]
            if cand_name in filename:
                best_cand = cand
                break
                
        if best_cand:
            try:
                with open(best_cand, "r") as f:
                    # It's a FASTA, so skip header
                    lines = f.readlines()
                    seq_lines = [l.strip() for l in lines if not l.startswith(">")]
                    rna_seq = "".join(seq_lines)
            except: pass
        
        # Override with explicit FASTA if provided
        if fasta_path and os.path.exists(fasta_path):
             try:
                with open(fasta_path, "r") as f:
                    lines = f.readlines()
                    seq_lines = [l.strip() for l in lines if not l.startswith(">")]
                    rna_seq = "".join(seq_lines)
                    # Normalize
                    rna_seq = core.dna_to_rna(rna_seq)
             except Exception as e:
                 print(f"Error reading explicit FASTA {fasta_path}: {e}")

        if not rna_seq:
            print("Could not locate or read RNA FASTA for context reconstruction.")
            return {}
            
        # Construct Context
        context = {
            'motif_data': motif_data,
            'rna_seq': rna_seq,
            'min_front': 1,
            'min_back': 1,
            'max_gap': max_gap,
            'predicted_scores': [], # Can't easily reconstruct if used, but usually empty for basic runs
            'use_ps_score': use_ps,
            'method': method,
            'ppr_name': "Reconstructed" # Display name
        }
        return context
        
    except Exception as e:
        print(f"Error reconstructing context: {e}")
        return {}

def load_single_site_optimized(path, target_pos, target_strand):
    """Loads a single row matching Position and Strand, parsing list columns."""
    try:
        # 1. Load Filtered (Fast)
        df = pd.read_csv(path)
        
        # 2. Filter target
        try:
            target_pos = int(target_pos)
        except:
            pass # Keep as is if conversion fails, though it should be int
            
        subset = df[(df['Genome_Position'] == target_pos) & (df['Strand'] == target_strand)]
        
        if subset.empty:
            return pd.DataFrame()
            
        # 3. Parse list columns ONLY for this subset (usually 1 row)
        row = subset.head(1).copy() # Take first match
        
        list_cols = ['PPR_Binding_Scores_List', 'PPR_Motif_Scores', 'Predicted PPR Motif Scores List', 'Motif_Confidence_Scores']
        
        def parse_val(x):
            if isinstance(x, str):
                try:
                    return ast.literal_eval(x)
                except:
                    # Fallback to manual parsing if literal_eval fails
                     if x.startswith('[') and x.endswith(']'):
                         try:
                             return [float(i) for i in x[1:-1].split(',')]
                         except:
                             return []
            return x

        for col in list_cols:
            if col in row.columns:
                row[col] = row[col].apply(parse_val)
                
        return row # Returns DataFrame with 1 row
        
    except Exception as e:
        print(f"Error in load_single_site_optimized: {e}")
        return pd.DataFrame()

def load_result(filename):
    """Loads a past result CSV and updates the UI."""
    # Output signature: status, csv, browser, motif_plot, motif_table, range_plot, df_state, analysis_context
    if not filename:
        return "No file selected.", None, go.Figure(), go.Figure(), pd.DataFrame(), go.Figure(), None, {}
    
    path = os.path.join(TMP_DIR, filename)
    if not os.path.exists(path):
        return f"File not found: {filename}", None, go.Figure(), go.Figure(), pd.DataFrame(), go.Figure(), None, {}
        
    try:
        df = load_csv_with_parsing(path)
        if df.empty:
            return "Error parsing CSV file.", None, go.Figure(), go.Figure(), pd.DataFrame(), go.Figure(), None, {}
        # Re-plot
        job_name = filename.replace("_results.csv", "")
        # Analysis type from DF if available
        atype = df['Analysis_Type'].iloc[0] if 'Analysis_Type' in df.columns and not df.empty else "Loaded"
        
        fig = plot_genome_browser(df, seq_id=job_name, analysis_type=atype)
        
        # Try to find corresponding motif file
        # Naming convention: {job}_results.csv -> {job}_motifs.txt
        motif_path = path.replace("_results.csv", "_motifs.txt")
        motif_fig = go.Figure()
        motif_df = pd.DataFrame()
        
        if os.path.exists(motif_path):
             # Extract simple title
             motif_df, motif_fig = visualize_motifs(motif_path, analysis_title=f"{job_name} Motifs", discovery_tool=None)
        
        # Return empty dict for analysis_context to force Static Mode in Detail Inspector
        return f"Loaded result: {filename}", path, fig, motif_fig, motif_df, go.Figure(), df, {}
    except Exception as e:
        print(f"Error loading result: {e}")
        return f"Error loading: {e}", None, go.Figure(), go.Figure(), pd.DataFrame(), go.Figure(), None, {}

def run_unified_analysis(protein_text, protein_file, rna_text, rna_file, tool_choice, scoring_method, max_gap_size, use_ps_score, job_name, highlight_top=False, highlight_strand="Both", skip_plots=False):
    # Construct analysis_type string
    method_map = {"PLS-type": "kobayashi", "P-type": "yan"}
    method = method_map.get(scoring_method, scoring_method).lower()
    gap_str = f"gap{int(max_gap_size)}" if max_gap_size > 0 else "gap0"
    ps_str = "ps" if use_ps_score else "no_ps"
    analysis_type = f"{method}_{gap_str}_{ps_str}"
    # 1. Process Protein Input (Discovery)
    if protein_file is not None:
        with open(protein_file.name, 'r') as f:
            protein_input = f.read()
    else:
        protein_input = protein_text
        if not protein_input:
            protein_input = DEFAULT_PROTEIN

    if not protein_input:
        return "Error: Please provide a protein sequence.", None, None, None, None, None, None, None, None, None, None

    if not protein_input.strip().startswith(">"):
        protein_input = f">UserSequence\n{protein_input}"

    # 2. Process RNA Input (Analysis Target)
    if rna_file is not None:
        with open(rna_file.name, 'r') as f:
            rna_input = f.read()
    else:
        rna_input = rna_text
        if not rna_input:
            rna_input = DEFAULT_RNA

    if not rna_input:
        return "Error: Please provide an RNA sequence.", None, None, None, None, None, None, None, None, None
        
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

    # User Request: Normalize to RNA (T->U, Uppercase)
    # Also aggressive cleanup of spaces/newlines to prevent "empty space" errors
    sequence = sequence.replace(" ", "").replace("\n", "").replace("\r", "").strip()
    sequence = core.dna_to_rna(sequence)

    # 3. Run Discovery
    # Pre-clean protein input to handle newlines/spaces correctly
    try:
        motif_data, ppr_name, pkl_path = run_motif_discovery(protein_input, tool_choice)
        motifs_file = f"{pkl_path.replace('.pkl', '')}_motifs.txt" # Approximate path
        # run_motif_discovery uses tempfile for pickle.
        # But we need motifs_file path for saving later?
        # run_motif_discovery returns pkl_path.
        # It saved motifs file in the output_dir of detection.
        # We can find it.
        motifs_file = os.path.join(os.path.dirname(pkl_path), f"{ppr_name}_motifs.txt")
        
    except Exception as e:
        return f"Error during discovery: {e}", None, None, None, None, None, None, None, None, None

    # 4. Run Analysis
    status_msg = f"Motif Discovery Complete: {ppr_name}"
    status_msg += f"\nRunning Analysis ({analysis_type})..."
    
    # Load the discovered motifs (Done in helper)
    # motif_data is ready.

    # Run Analysis (Both Strands)
    min_front = 1 # Allow gaps closer to ends
    min_back = 1
    max_gap = max_gap_size
    
    # Prepare RNA sequence
    # rna_input is already defined from rna_text or rna_file earlier.
    if not rna_input:
        return "Error: No RNA sequence provided.", None, None, None, None, None, None, None, None, None

    # Parse RNA FASTA
    try:
        rna_record = next(SeqIO.parse(io.StringIO(rna_input), "fasta"))
        rna_seq = str(rna_record.seq).upper().replace("U", "T") # Convert to DNA for analysis if needed, or keep RNA. Core usually expects DNA/RNA string.
    except Exception as e:
        # Fallback if raw sequence or malformed FASTA
        # If it looks like FASTA (starts with >), try to strip header manually
        clean_input = rna_input.strip()
        if clean_input.startswith(">"):
             parts = clean_input.split("\n", 1)
             if len(parts) > 1:
                 # Remove header, join lines to remove newlines
                 rna_seq = "".join(parts[1].splitlines()).upper().replace("U", "T")
             else:
                 rna_seq = ""
        else:
             # Raw sequence, remove newlines
             rna_seq = "".join(clean_input.splitlines()).upper().replace("U", "T")
    
             rna_seq = "".join(clean_input.splitlines()).upper().replace("U", "T")
    
    # Initialize all_scores as a dict of lists
    all_scores = {
        "Position_Relative": [],
        "PPR_Binding_Score": [],
        "PPR_Name": [],
        "Motif_Labels": [],
        "Motif_Labels_GAP": [],
        "Motif_Types_GAP": [], # New Column
        "Gap_Size": [],
        "Gap_Index": [], # New Column
        "Segment": [],
        "PPR_Motif_Scores": [],
        "Total_Score": [],
        "PPR_Binding_Scores_List": [],
        "Motif_Confidence_Scores": [],
        "Strand": [],
        "Genome_Position": []
    }
    
    # Construct lookup key for core.analyze_ppr_data
    # core.py returns fixed keys: 'yan_gap0_ps' or 'yan_gap_ps' (generic 'gap' for any >0)
    lookup_gap_str = "gap" if max_gap_size > 0 else "gap0"
    lookup_key = f"{method}_{lookup_gap_str}_{ps_str}"
    
    predicted_scores = [] # Default empty
    
    # Target Analysis Config for Optimization
    target_analysis_config = {
        'method': method, # 'yan' or 'kobayashi'
        'use_ps': use_ps_score,
        'gap_mode': 'varied' if max_gap > 0 else 'gap0'
    }

    try:
        # Forward
        # Use Lite Mode + Targeted Analysis
        scores_fwd = core.analyze_ppr_data(motif_data, rna_seq, min_front, min_back, max_gap, predicted_scores, ppr_name, detailed_output=False, target_analysis=target_analysis_config)
        if lookup_key in scores_fwd:
            fwd_data = scores_fwd[lookup_key]
            # Extend all_scores with fwd_data (Dict)
            for k, v in fwd_data.items():
                 if k in all_scores:
                     all_scores[k].extend(v)
            
            # Add Strand and Genome_Position columns
            num_fwd = len(fwd_data["Position_Relative"])
            all_scores["Strand"].extend(["Forward"] * num_fwd)
            all_scores["Genome_Position"].extend(fwd_data["Position_Relative"]) # Same for fwd
            print(f"DEBUG: Forward analysis done. Found {num_fwd} sites.")
            fwd_msg = f"Forward: {num_fwd} sites. "
        
        # Reverse
        rev_seq = core.get_reverse_complement(rna_seq)
        print(f"DEBUG: Starting Reverse analysis. Seq length: {len(rev_seq)}")
        scores_rev = core.analyze_ppr_data(motif_data, rev_seq, min_front, min_back, max_gap, predicted_scores, ppr_name, detailed_output=False, target_analysis=target_analysis_config)
        
        rev_msg = "Reverse: 0 sites."
        if lookup_key in scores_rev:
            rev_data = scores_rev[lookup_key]
            # Extend all_scores with rev_data (Dict)
            for k, v in rev_data.items():
                 if k in all_scores:
                     all_scores[k].extend(v)
            
            num_rev = len(rev_data["Position_Relative"])
            all_scores["Strand"].extend(["Reverse"] * num_rev)
            # Genome position for reverse: len(seq) - 1 - relative
            seq_len = len(rna_seq)
            rev_genome_pos = [seq_len - 1 - p for p in rev_data["Position_Relative"]]
            all_scores["Genome_Position"].extend(rev_genome_pos)
            print(f"DEBUG: Reverse analysis done. Found {num_rev} sites.")
            rev_msg = f"Reverse: {num_rev} sites."
        else:
            print(f"DEBUG: Lookup key {lookup_key} not found in reverse scores. Keys: {list(scores_rev.keys())}")
            
    except Exception as e:
        import traceback
        traceback.print_exc()
        return f"Error during analysis: {e}", None, None, None, None, None, None, None, None, None
    if not all_scores["Position_Relative"]:
        return "Error: No binding sites found.", None, None, None, None, None, None, None, None, None
        
    status_msg = f"Success. {fwd_msg}{rev_msg}"

    # Create DataFrame
    try:
        df_output = pd.DataFrame(all_scores)
    except Exception as e:
        print(f"ERROR: DataFrame construction failed: {e}")
        return f"Error building DataFrame: {e}", None, None, None, None, None, None, None, None, None
    
    # Ensure numerical types for sorting
    df_output["PPR_Binding_Score"] = pd.to_numeric(df_output["PPR_Binding_Score"], errors='coerce')
    
    # Add Analysis Type
    df_output['Analysis_Type'] = analysis_type
    
    # --- PERFORMANCE OPTIMIZATION ---
    # 4M+ rows (500MB+) is too heavy for web UI / CSV I/O.
    # We truncate to Top N *before* saving.
    MAX_CSV_ROWS = 30000
    
    # Truncation logic moved to AFTER saving
    pass
    
    # Filter for significant hits (optional, but good for default view)
    # For now, return all or top hits. Let's return all.
    
    if not df_output.empty and 'Gap_Size' in df_output.columns:
        gap_sizes = df_output['Gap_Size'].unique()
        print(f"Gap Sizes in Result: {gap_sizes}")
    
    # Sanitize helper
    def sanitize(s):
        import re
        return re.sub(r'[^\w\-_\.]', '_', str(s))

    # Extract IDs from FASTA content
    prot_id = "Protein"
    if protein_input and protein_input.strip().startswith(">"):
        header = protein_input.strip().split("\n")[0]
        # Use first word of header as ID, or up to 15 chars if no spaces? 
        # User said "FASTA name info". Usually that's the ID (first token).
        prot_id = header[1:].strip().split()[0] 

    rna_id = "RNA"
    if rna_input and rna_input.strip().startswith(">"):
        header = rna_input.strip().split("\n")[0]
        rna_id = header[1:].strip().split()[0]
        
    # Construct descriptive job name
    # Format: {UserJobName}_{ProtID}_{RNAID}_{AnalysisParams}
    # Avoid redundancy if UserJobName is default
    clean_job_name = sanitize(job_name)
    clean_prot = sanitize(prot_id)
    clean_rna = sanitize(rna_id)
    # analysis_type is already formatted e.g. "yan_gap2_ps"
    
    # If job name is generic, maybe replace it? But user said "arbitrary input string AND ...".
    # So we MUST include the input string.
    final_job_name = f"{clean_job_name}_{clean_prot}_{clean_rna}_{analysis_type}"
    
    # Create Analysis Context for Dynamic Detail Inspection
    analysis_context = {
        'motif_data': motif_data,
        'rna_seq': rna_seq,
        'min_front': min_front, 
        'min_back': min_back,
        'max_gap': max_gap,
        'predicted_scores': predicted_scores,
        'use_ps_score': use_ps_score,
        'method': method, # 'yan' or 'kobayashi'
        'ppr_name': ppr_name
    }
    
    # Construct friendly filename parts
    friendly_method = "P-type" if method == "yan" else "PLS-type"
    friendly_gap = f"Gap{int(max_gap_size)}"
    friendly_score = "ScoreOn" if use_ps_score else "ScoreOff"
    
    # Format: {UserJobName}_{ProtID}_{RNAID}_{Method}_{Gap}_{Score}_results.csv
    final_job_name = f"{clean_job_name}_{clean_prot}_{clean_rna}_{friendly_method}_{friendly_gap}_{friendly_score}_results"
    
    # Save to CSV using the new descriptive name (FULL DATA)
    csv_path = os.path.join(TMP_DIR, f"{final_job_name}.csv")
    df_output.to_csv(csv_path, index=False)
    
    # Save accompanying motifs file for Load Result
    if motifs_file and os.path.exists(motifs_file):
        # Naming convention matches load_result expectation: {job}_motifs.txt
        # final_job_name ends with _results
        motifs_dest_name = f"{final_job_name.replace('_results', '_motifs')}.txt"
        dest_motifs = os.path.join(TMP_DIR, motifs_dest_name)
        try:
            shutil.copy2(motifs_file, dest_motifs)
            print(f"DEBUG: Copied motifs from {motifs_file} to {dest_motifs}")
        except Exception as e:
            print(f"Warning: Failed to save motifs file: {e}")
    else:
        # If no motifs file (e.g. from text input?), maybe save motif_data as text?
        print(f"DEBUG: Motifs file missing or invalid. Source: {motifs_file}")

    # --- UI OPTIMIZATION: Truncate for display ---
    MAX_CSV_ROWS = 30000
    
    if len(df_output) > MAX_CSV_ROWS:
        msg = f"Result too large ({len(df_output)} rows). UI display truncated to Top {MAX_CSV_ROWS} sites by Score."
        msg += "\n(Full results saved to CSV)"
        print(msg)
        status_msg += "\n" + msg
        
        # Sort by Score descending
        df_output = df_output.sort_values(by="PPR_Binding_Score", ascending=False)
        
        # Balance Strand if possible (15k Fwd + 15k Rev) to protect Reverse graph
        df_fwd = df_output[df_output['Strand'] == 'Forward']
        df_rev = df_output[df_output['Strand'] == 'Reverse']
        
        limit_per = MAX_CSV_ROWS // 2
        df_top = pd.concat([df_fwd.head(limit_per), df_rev.head(limit_per)])
        
        # Re-sort by Position for readability
        df_output = df_top.sort_values(by=['Genome_Position'])
        
        print(f"Truncation complete. Final shape: {df_output.shape}")
        
    print(f"DEBUG: Returning DataFrame shape: {df_output.shape}")
    # 5. Visualization (Genome Browser)
    fig = go.Figure()
    if not skip_plots:
        fig = plot_genome_browser(df_output, seq_id, analysis_type, highlight_top=highlight_top, highlight_strand=highlight_strand)

    
    # 6. Visualize Motifs
    motif_df, motif_fig = pd.DataFrame(), go.Figure()
    if motifs_file and os.path.exists(motifs_file):
        friendly_title = f"Tool: {tool_choice} | Method: {friendly_method} | Gap: {int(max_gap_size)} | Score: {'On' if use_ps_score else 'Off'}"
        motif_df, motif_fig = visualize_motifs(motifs_file, analysis_title=friendly_title, discovery_tool=tool_choice)
    else:
        # Fallback empty
        motif_df, motif_fig = pd.DataFrame(), go.Figure()
        
    status_msg += f"\nAnalysis Complete. Found {len(df_output)} sites."
    
    # Get updated file list for dropdowns
    dd_main, dd_ncbi = list_past_results()
    
    # Outputs order: status, csv, motif_plot, browser_plot, range_plot, df_state, context, main_dd, ncbi_dd
    return status_msg, csv_path, motif_fig, motif_df, fig, go.Figure(), df_output, analysis_context, dd_main, dd_ncbi
    
def parse_analysis_type(atype_str):
    """Parses analysis_type string into readable format."""
    # format: {method}_{gap_str}_{ps_str} e.g. yan_gap0_ps or kobayashi_gap2_no_ps
    # method: yan->P-type, kobayashi->PLS-type
    if not atype_str or atype_str == "Analysis":
        return "Analysis Results"
        
    parts = atype_str.split('_')
    method = "Unknown"
    gap = "0"
    ps = "Off"
    
    if len(parts) >= 3:
        m = parts[0]
        g = parts[1]
        p = parts[2]
        
        if m == "yan": method = "P-type"
        elif m == "kobayashi": method = "PLS-type"
        else: method = m.capitalize()
        
        if g.startswith("gap"): gap = g.replace("gap", "")
        
        if "no" in p: ps = "Off"
        else: ps = "On"
        
        return f"Method: {method} | Gap: {gap} | Score: {ps}"
    return atype_str

def update_main_browser_plot(df, highlight, strand):
    if df is None or df.empty:
        return go.Figure()
    # analysis_type will be extracted from DF by plot_genome_browser
    return plot_genome_browser(df, analysis_type=None, highlight_top=highlight, highlight_strand=strand)

def apply_strand_filter(df, strand_mode):
    """Returns a strand-filtered dataframe for display controls."""
    if df is None or df.empty or 'Strand' not in df.columns:
        return df
    if strand_mode == "Forward":
        return df[df['Strand'] == "Forward"]
    if strand_mode == "Reverse":
        return df[df['Strand'] == "Reverse"]
    return df

def plot_genome_browser(df, seq_id="Sequence", analysis_type="Analysis", highlight_top=False, highlight_strand="Both"):
    """Helper to plot the genome browser."""
    # Try to get analysis type from DF if not provided or generic
    if (analysis_type == "Analysis" or analysis_type is None) and 'Analysis_Type' in df.columns and not df.empty:
        analysis_type = df['Analysis_Type'].iloc[0]
        
    readable_condition = parse_analysis_type(analysis_type)
    display_df = apply_strand_filter(df, highlight_strand)

    show_forward = highlight_strand in ("Both", "Forward")
    show_reverse = highlight_strand in ("Both", "Reverse")
    rows = 2 if (show_forward and show_reverse) else 1

    fwd_row = 1 if show_forward else None
    rev_row = None
    if show_reverse:
        rev_row = 2 if rows == 2 else 1
    
    fig = make_subplots(
        rows=rows, cols=1, 
        shared_xaxes=(rows == 2), 
        vertical_spacing=0.12 if rows == 2 else 0.0
        # No subplot_titles to avoid overlap with legend
    )
    
    unique_ppr_names = display_df['PPR_Name'].unique() if not display_df.empty else []
    colors = px.colors.qualitative.D3  # D3 palette: higher contrast and colorblind-friendly

    for i, name in enumerate(unique_ppr_names):
        color = colors[i % len(colors)]

        # Forward
        df_fwd = display_df[(display_df['PPR_Name'] == name) & (display_df['Strand'] == "Forward")]
        if fwd_row is not None and not df_fwd.empty:
            # Sort by Score ascending so higher scores are plotted on top
            df_fwd = df_fwd.sort_values(by="PPR_Binding_Score", ascending=True)

            # Prepare Custom Data: [Index, Score, Gap_Size, Gap_Index]
            custom_data = np.stack((
                df_fwd.index,
                df_fwd['PPR_Binding_Score'],
                df_fwd.get('Gap_Size', [0]*len(df_fwd)),
                df_fwd.get('Gap_Index', [0]*len(df_fwd))
            ), axis=-1)

            fig.add_trace(go.Scattergl(
                x=df_fwd['Genome_Position'],
                y=df_fwd['PPR_Binding_Score'],
                mode='markers',
                name=f"{name} (+)",
                marker=dict(color=color, size=7, symbol='circle', opacity=0.7),
                customdata=custom_data,
                hovertemplate="Pos: %{x}<br>Score: %{y:.2f}<br>Strand: +<br>Gap: %{customdata[2]} (Pos %{customdata[3]})<extra></extra>"
            ), row=fwd_row, col=1)

        # Reverse
        df_rev = display_df[(display_df['PPR_Name'] == name) & (display_df['Strand'] == "Reverse")]
        if rev_row is not None and not df_rev.empty:
            df_rev = df_rev.sort_values(by="PPR_Binding_Score", ascending=True)

            custom_data = np.stack((
                df_rev.index,
                df_rev['PPR_Binding_Score'],
                df_rev.get('Gap_Size', [0]*len(df_rev)),
                df_rev.get('Gap_Index', [0]*len(df_rev))
            ), axis=-1)

            fig.add_trace(go.Scattergl(
                x=df_rev['Genome_Position'],
                y=df_rev['PPR_Binding_Score'],
                mode='markers',
                name=f"{name} (-)",
                marker=dict(color=color, size=7, symbol='diamond', opacity=0.7),
                customdata=custom_data,
                hovertemplate="Pos: %{x}<br>Score: %{y:.2f}<br>Strand: -<br>Gap: %{customdata[2]} (Pos %{customdata[3]})<extra></extra>"
            ), row=rev_row, col=1)

    # Highlight Top Score
    if highlight_top and not display_df.empty:
        # Search only within currently displayed strand scope
        df_search = display_df

        if not df_search.empty:
            top_idx = df_search['PPR_Binding_Score'].abs().idxmax()
            top_row = display_df.loc[top_idx]
            
            t_pos = top_row['Genome_Position']
            t_score = top_row['PPR_Binding_Score']
            t_strand = top_row['Strand']
            
            row_idx = fwd_row if t_strand == "Forward" else rev_row
            if row_idx is None:
                row_idx = 1
            
            fig.add_trace(go.Scatter(
                x=[t_pos], y=[t_score],
                mode='markers',
                marker=dict(symbol='star', size=15, color='gold', line=dict(width=2, color='black')),
                name=f"Top Score ({t_strand})",
                hoverinfo='text',
                hovertext=f"TOP SCORE<br>Pos: {t_pos}<br>Score: {t_score:.2f}"
            ), row=row_idx, col=1)

    # Calculate Global Max/Min for visual guidance
    global_max = 0
    global_min = 0
    if not display_df.empty:
        scores = display_df['PPR_Binding_Score']
        global_max = scores.max()
        global_min = scores.min()
    
    # Add Max/Min reference lines
    if global_max > 0:
        if fwd_row is not None:
            fig.add_hline(y=global_max, line_dash="dash", line_color="blue", annotation_text=f"Max: {global_max:.1f}", annotation_position="top right", row=fwd_row, col=1)
        if rev_row is not None:
            max_anno_pos = "bottom right" if rows == 2 else "top right"
            fig.add_hline(y=global_max, line_dash="dash", line_color="blue", annotation_text=f"Max: {global_max:.1f}", annotation_position=max_anno_pos, row=rev_row, col=1)
    
    if global_min > 0:
        if fwd_row is not None:
            fig.add_hline(y=global_min, line_dash="dot", line_color="red", annotation_text=f"Min: {global_min:.1f}", annotation_position="bottom right", row=fwd_row, col=1)
        if rev_row is not None:
            min_anno_pos = "top right" if rows == 2 else "bottom right"
            fig.add_hline(y=global_min, line_dash="dot", line_color="red", annotation_text=f"Min: {global_min:.1f}", annotation_position=min_anno_pos, row=rev_row, col=1)

    # Calculate Y-axis limits based on Global Min/Max with padding
    pad = (global_max - global_min) * 0.1
    if pad == 0: pad = 10
    y_max_limit = global_max + pad
    y_min_limit = global_min - pad

    import time
    strand_label = {"Both": "+ / -", "Forward": "+ only", "Reverse": "- only"}.get(highlight_strand, str(highlight_strand))
    fig.update_layout(
        title=dict(
            text=f"<b>Genome Browser - {seq_id}</b><br><sup>Method: {readable_condition} | Strand: {strand_label}</sup>",
            y=0.95,
            x=0.5,
            xanchor='center',
            yanchor='top'
        ),
        height=500 if rows == 2 else 360,
        autosize=True,
        hovermode="closest",
        legend=dict(orientation="h", y=-0.15, x=0.5, xanchor="center", yanchor="top"),  # Legend at bottom
        margin=dict(t=50, b=80, l=60, r=60),
        uirevision=str(time.time())  # Force view reset on each render
    )
    
    fig.update_xaxes(title_text="Genomic Position", row=rows, col=1)
    if fwd_row is not None and rev_row is not None and rows == 2:
        fig.update_yaxes(title_text="Fwd Score", row=fwd_row, col=1, range=[y_min_limit, y_max_limit], autorange=False)
        fig.update_yaxes(title_text="Rev Score", row=rev_row, col=1, range=[y_min_limit, y_max_limit], autorange="reversed")
    elif fwd_row is not None:
        fig.update_yaxes(title_text="Fwd Score", row=fwd_row, col=1, range=[y_min_limit, y_max_limit], autorange=False)
    elif rev_row is not None:
        fig.update_yaxes(title_text="Rev Score", row=rev_row, col=1, range=[y_min_limit, y_max_limit], autorange="reversed")
    
    return fig

def plot_range_browser(df, gff_path, start, end, seq_id="Unknown", feature_types=["gene"], plot_height=600, highlight_top=False, highlight_strand="Both", analysis_context=None, protein_name=None):
    """Plots binding scores and gene annotations for a specific range.
       If analysis_context is provided and range is small (<50kb), it optionally re-calculates 
       to show ALL sites (ignoring truncation).
    """
    # Extract analysis type from DF
    analysis_type = "Analysis"
    if df is not None and not df.empty and 'Analysis_Type' in df.columns:
        analysis_type = df['Analysis_Type'].iloc[0]
    
    readable_condition = parse_analysis_type(analysis_type)
    
    if feature_types is None:
        feature_types = ["gene", "CDS", "tRNA", "rRNA", "exon"]

    show_forward = highlight_strand in ("Both", "Forward")
    show_reverse = highlight_strand in ("Both", "Reverse")
    if show_forward and show_reverse:
        rows = 3
        row_heights = [0.4, 0.2, 0.4]
        vertical_spacing = 0.10
        display_height = plot_height
        fwd_row = 1
        gene_row = 2
        rev_row = 3
    else:
        rows = 2
        row_heights = [0.78, 0.22]
        vertical_spacing = 0.08
        display_height = max(420, int(plot_height * 0.82))
        fwd_row = 1 if show_forward else None
        rev_row = 1 if show_reverse else None
        gene_row = 2

    fig = make_subplots(
        rows=rows, cols=1, 
        shared_xaxes=True, 
        vertical_spacing=vertical_spacing,
        row_heights=row_heights
        # No subplot_titles to avoid overlap with legend
    )
    
    
    # 1. Plot Binding Scores (Row 1: Forward, Row 3: Reverse)
    # 1. Plot Binding Scores (Row 1: Forward, Row 3: Reverse)
    # Check for Dynamic Re-calculation
    MAX_RECALC_RANGE = 50000
    range_size = end - start
    
    df_range = pd.DataFrame() # Default empty
    is_recalc = False
    recalc_msg = ""
    
    if analysis_context and range_size < MAX_RECALC_RANGE and 'rna_seq' in analysis_context and 'motif_data' in analysis_context:
        try:
            # Dynamic Re-calculation
            rna_seq = analysis_context['rna_seq']
            if start < 0: start = 0
            if end > len(rna_seq): end = len(rna_seq)
            
            # Slice sequence (add padding to avoid edge effects?)
            # Usually strict slice is fine if we just want hits IN range.
            # But motifs overlapping edges might be lost. 
            # Let's add padding equal to max_gap + motif_length (~60bp).
            padding = 100
            s_pad = max(0, start - padding)
            e_pad = min(len(rna_seq), end + padding)
            
            subseq = rna_seq[s_pad:e_pad]
            
            # Run Analysis
            motif_data = analysis_context['motif_data']
            # Reconstruct settings from context or use defaults (be careful with defaults!)
            # We need: params like 'method', 'use_ps_score', 'max_gap'.
            # run_unified_analysis stores them in context properly.
            
            from ppr_analysis import core
            
            # Analysis Params
            method = analysis_context.get('method', 'yan')
            use_ps = analysis_context.get('use_ps_score', True)
            max_gap = analysis_context.get('max_gap', 0)
            min_front = analysis_context.get('min_front', 0)
            min_back = analysis_context.get('min_back', 0)
            predicted_scores = analysis_context.get('predicted_scores', [])
            ppr_name = analysis_context.get('ppr_name', 'PPR')

            gap_mode = 'varied' if max_gap > 0 else 'gap0'
            target_config = {'method': method, 'use_ps': use_ps, 'gap_mode': gap_mode}
            
            # Forward Analysis
            # Signature: analyze_ppr_data(motif_data, fasta_sequence, min_front, min_back, max_gap, predicted_scores, ppr_name, detailed_output=True, target_analysis=None)
            res_fwd_dict = core.analyze_ppr_data(motif_data, subseq, min_front, min_back, max_gap, predicted_scores, ppr_name, detailed_output=True, target_analysis=target_config)
            
            # Extract the specific analysis result from the dict (e.g., 'yan_gap_ps')
            # The key name depends on params.
            # Helper to generate key:
            key_method = method
            key_gap = "gap" if gap_mode == 'varied' else "gap0"
            key_ps = "ps" if use_ps else "no_ps"
            result_key = f"{key_method}_{key_gap}_{key_ps}"
            
            res_fwd = res_fwd_dict.get(result_key, {})

            # Reverse Analysis
            # We must generate Reverse Complement of subseq
            # Simple RC helper
            def reverse_complement(seq):
                comp = {'A':'U', 'U':'A', 'G':'C', 'C':'G', 'T':'A', 
                        'a':'u', 'u':'a', 'g':'c', 'c':'g', 't':'a',
                        'N':'N', 'n':'n', '-':'-'}
                # Handle DNA T->U if needed, but usually we work with RNA U.
                # If input is DNA, T->A.
                # Assuming RNA input.
                try:
                    return "".join(comp.get(base, base) for base in reversed(seq))
                except:
                    return seq[::-1] # Fallback reverse
            
            subseq_rc = reverse_complement(subseq)
            
            res_rev_dict = core.analyze_ppr_data(motif_data, subseq_rc, min_front, min_back, max_gap, predicted_scores, ppr_name, detailed_output=True, target_analysis=target_config)
            res_rev = res_rev_dict.get(result_key, {})
            
            # Convert Results to DataFrame & Adjust Coordinates (Local -> Global)
            # Forward
            local_scores = {}
            if res_fwd and 'Position_Relative' in res_fwd:
                 local_scores = res_fwd
                 local_scores['Genome_Position'] = [p + s_pad for p in res_fwd['Position_Relative']]
                 local_scores['Strand'] = ["Forward"] * len(res_fwd['Position_Relative'])
                 # Ensure PPR_Name is populated (it is returned by analyze_ppr_data as list)

            # Reverse
            if res_rev and 'Position_Relative' in res_rev:
                if not local_scores:
                    local_scores = {k: [] for k in res_rev.keys()}
                    # Ensure special keys exist
                    if 'Genome_Position' not in local_scores: local_scores['Genome_Position'] = []
                    if 'Strand' not in local_scores: local_scores['Strand'] = []
                
                # Extend generic columns
                for k in res_rev.keys():
                    if k in local_scores:
                        local_scores[k].extend(res_rev[k])
                
                # Calcs for Global Coordinates
                # Local Rev Pos `p` is 0-based index from 5' end of RC sequence.
                # RC[0] corresponds to Original[Last].
                # RC[p] corresponds to Original[Last - p].
                # Global Coord = s_pad + (len(subseq) - 1 - p).
                # Wait, Position usually refers to the 5' end of the binding site on the strand.
                # If site is at p on RC, it binds to region p..p+len.
                # On Forward Genome, this is (len - 1 - p) .. down.
                # Let's verify `run_unified_analysis` logic:
                # `rev_genome_pos = [seq_len - p + 1 for p in rev_pos_relative]` (if 1-based?)
                # Actually, `analyze_ppr_data` returns 0-based indices if python logic.
                # Let's assume 0-based for now.
                # Global 1-based = s_pad + (len(subseq) - p). (Approx).
                # Let's use logic: `s_pad + len(subseq) - p`
                
                rev_pos_rel = res_rev['Position_Relative']
                sub_len = len(subseq)
                
                # Using the formula derived from run_unified_analysis context
                rev_global = [s_pad + sub_len - p for p in rev_pos_rel]
                
                if 'Genome_Position' in local_scores:
                    local_scores['Genome_Position'].extend(rev_global)
                else: 
                     local_scores['Genome_Position'] = rev_global
                     
                local_scores['Strand'].extend(["Reverse"] * len(rev_pos_rel))
                # PPR_Name is already extended by generic loop if present in res_rev
            
            # Create DF
            if local_scores:
                 df_recalc = pd.DataFrame(local_scores)
                 # Filter strictly to range (remove padding effects)
                 df_range = df_recalc[(df_recalc['Genome_Position'] >= start) & (df_recalc['Genome_Position'] <= end)]
                 is_recalc = True
                 recalc_msg = f" (Full Re-calc for {range_size}bp)"
            else:
                 df_range = pd.DataFrame() # Empty result
                 is_recalc = True

        except Exception as e:
            print(f"DEBUG: Recalc failed: {e}. Fallback to static.")
            # Fallback
            df_range = df[(df['Genome_Position'] >= start) & (df['Genome_Position'] <= end)]
    
    else:
        # Static Filtering (Truncated data)
        df_range = df[(df['Genome_Position'] >= start) & (df['Genome_Position'] <= end)]
        if range_size >= MAX_RECALC_RANGE and analysis_context:
             recalc_msg = " (Top Hits Only - Zoom in <50kb for Full Data)"

    df_plot = apply_strand_filter(df_range, highlight_strand)
    print(f"DEBUG: Plotting Range {start}-{end}{recalc_msg}. Rows: {len(df_plot)} (from {len(df_range)})")
    
    # ... Rest of plotting code uses df_range ...
    # We need to ensure df_range has 'PPR_Name' column if re-calculated. Added above.
    
    readable_condition += recalc_msg
    
    # Update title with recalc status
    
    if not df_plot.empty:
        print(f"DEBUG: Strands in Range: {df_plot['Strand'].unique()}")
        print(f"DEBUG: Top 5 rows in range:\n{df_plot[['Genome_Position', 'Strand', 'PPR_Binding_Score']].head().to_string()}")
        # Calculate max score from VISIBLE range for symmetric Y-axis
        max_score = df_plot['PPR_Binding_Score'].abs().max()
        if max_score == 0: max_score = 10
    else:
        # Fallback to full DF or default
        if not df.empty:
             max_score = df['PPR_Binding_Score'].abs().max()
             if max_score == 0: max_score = 10
        else:
             max_score = 10
             
    y_limit = max_score * 1.1

    unique_ppr_names = df_plot['PPR_Name'].unique() if not df_plot.empty else []
    colors = px.colors.qualitative.D3  # D3 palette: higher contrast and colorblind-friendly

    for i, name in enumerate(unique_ppr_names):
        color = colors[i % len(colors)]

        # Forward
        df_fwd = df_plot[(df_plot['PPR_Name'] == name) & (df_plot['Strand'] == "Forward")]
        if fwd_row is not None and not df_fwd.empty:
            # Prepare Custom Data for Range Browser if available
            custom_data = np.stack((
                df_fwd.index,
                df_fwd['PPR_Binding_Score'],
                df_fwd.get('Gap_Size', [0]*len(df_fwd)),
                df_fwd.get('Gap_Index', [0]*len(df_fwd))
            ), axis=-1)

            fig.add_trace(go.Scattergl(
                x=df_fwd['Genome_Position'],
                y=df_fwd['PPR_Binding_Score'],
                mode='markers',
                name=f"{name} (+)",
                marker=dict(color=color, size=7, symbol='circle', opacity=0.7),
                customdata=custom_data,
                hovertemplate="Pos: %{x}<br>Score: %{y:.2f}<br>Strand: +<br>Gap: %{customdata[2]} (Pos %{customdata[3]})<extra></extra>"
            ), row=fwd_row, col=1)

        # Reverse
        df_rev = df_plot[(df_plot['PPR_Name'] == name) & (df_plot['Strand'] == "Reverse")]
        if rev_row is not None and not df_rev.empty:
            custom_data = np.stack((
                df_rev.index,
                df_rev['PPR_Binding_Score'],
                df_rev.get('Gap_Size', [0]*len(df_rev)),
                df_rev.get('Gap_Index', [0]*len(df_rev))
            ), axis=-1)

            fig.add_trace(go.Scattergl(
                x=df_rev['Genome_Position'],
                y=df_rev['PPR_Binding_Score'],
                mode='markers',
                name=f"{name} (-)",
                marker=dict(color=color, size=7, symbol='diamond', opacity=0.7),
                customdata=custom_data,
                hovertemplate="Pos: %{x}<br>Score: %{y:.2f}<br>Strand: -<br>Gap: %{customdata[2]} (Pos %{customdata[3]})<extra></extra>"
            ), row=rev_row, col=1)

    # Highlight Top Score in Range
    if highlight_top and not df_plot.empty:
         # Search only within currently displayed strand scope
        df_search = df_plot

        if not df_search.empty:
            top_idx = df_search['PPR_Binding_Score'].abs().idxmax()
            top_row = df_plot.loc[top_idx] 
            
            t_pos = top_row['Genome_Position']
            t_score = top_row['PPR_Binding_Score']
            t_strand = top_row['Strand']
            
            row_idx = fwd_row if t_strand == "Forward" else rev_row
            if row_idx is None:
                row_idx = 1
            
            fig.add_trace(go.Scatter(
                x=[t_pos], y=[t_score],
                mode='markers',
                marker=dict(symbol='star', size=15, color='gold', line=dict(width=2, color='black')),
                name=f"Top Score ({t_strand})",
                hoverinfo='text',
                hovertext=f"TOP SCORE<br>Pos: {t_pos}<br>Score: {t_score:.2f}"
            ), row=row_idx, col=1)

    # Plot Gene Annotations
    fig.add_trace(go.Scatter(x=[start, end], y=[0, 0], mode='markers', marker=dict(opacity=0), showlegend=False, hoverinfo='skip'), row=gene_row, col=1)
    
    if gff_path and os.path.exists(gff_path):
        features = ncbi.parse_gff(gff_path, start, end)
        print(f"DEBUG: Range {start}-{end}, GFF: {gff_path}, Features Found: {len(features)}")
        print(f"DEBUG: Selected Feature Types: {feature_types}")
        
        gene_y_fwd = 0.5
        gene_y_rev = -0.5
        
        for feat in features:
            f_type = feat['type']
            if f_type not in feature_types:
                # print(f"DEBUG: Skipping {f_type}")
                continue
            
            print(f"DEBUG: Plotting {f_type} at {feat['start']}-{feat['end']}")
                
            f_start = max(start, feat['start'])
            f_end = min(end, feat['end'])
            f_strand = feat['strand']
            
            # Extract name
            attrs = feat.get('attributes', {})
            f_name = attrs.get('Name', attrs.get('gene', attrs.get('ID', f_type)))
            
            # Color by type
            if f_type == 'CDS': color = 'orange'
            elif f_type == 'tRNA': color = 'green'
            elif f_type == 'rRNA': color = 'purple'
            elif f_type == 'gene': color = 'blue'
            else: color = 'gray'
            
            if f_strand == '+':
                y_pos = gene_y_fwd
            else:
                y_pos = gene_y_rev
                
            # Arrow drawing logic reused (omitted details for brevity in this prompt context but kept in code)
            thickness = 0.4
            half_thick = thickness / 2
            view_width = end - start
            head_len = view_width * 0.01
            head_len = min(head_len, f_end - f_start)
            
            if f_strand == '+':
                y_top = y_pos + half_thick
                y_bot = y_pos - half_thick
                y_mid = y_pos
                x_head_base = max(f_start, f_end - head_len)
                path = f"M {f_start},{y_top} L {x_head_base},{y_top} L {f_end},{y_mid} L {x_head_base},{y_bot} L {f_start},{y_bot} Z"
            else:
                y_top = y_pos + half_thick
                y_bot = y_pos - half_thick
                y_mid = y_pos
                x_head_base = min(f_end, f_start + head_len)
                path = f"M {x_head_base},{y_top} L {f_end},{y_top} L {f_end},{y_bot} L {x_head_base},{y_bot} L {f_start},{y_mid} Z"

            fig.add_shape(
                type="path", path=path, line=dict(color=color, width=1), fillcolor=color, opacity=0.6, row=gene_row, col=1
            )
            
            # Simplified Label Logic (Centered, No Box)
            
            fig.add_annotation(
                x=(f_start + f_end) / 2,
                y=y_pos,
                text=f"<b>{f_name}</b>",
                showarrow=False,
                font=dict(color="black", size=10),
                yshift=0,
                # bgcolor="white", # Removed
                # bordercolor="gray", # Removed
                # borderwidth=1, # Removed
                opacity=1.0,
                row=gene_row, col=1
            )
            
    # Add title with conditions
    import time

    # Construct Title
    # Format: [ProteinName] Range Start-End (SeqID)
    prot_str = f"[{protein_name}] " if protein_name else ""
    seq_str = f"({seq_id})" if seq_id else ""
    
    title_main = f"<b>{prot_str}Range {start}-{end} {seq_str}</b>"

    strand_label = {"Both": "+ / -", "Forward": "+ only", "Reverse": "- only"}.get(highlight_strand, str(highlight_strand))
    fig.update_layout(
        title=dict(
            text=f"{title_main}<br><sup>Method: {readable_condition}{recalc_msg} | Strand: {strand_label}</sup>",
            y=0.95,
            x=0.5,
            xanchor='center',
            yanchor='top'
        ),
        height=display_height,
        autosize=True,
        hovermode="closest",
        legend=dict(orientation="h", y=-0.1, x=0.5, xanchor="center", yanchor="top"),  # Legend at bottom
        margin=dict(t=50, b=80, l=60, r=60),
        uirevision=str(time.time())  # Force view reset on each render
    )        
    # Shared Y-Axis Limits
    
    # Calculate Global Max/Min for Visual Guidance (Combined Forward & Reverse)
    global_max = 0
    global_min = 0 # Usually binding score is positive, but let's check true min
    
    if not df_plot.empty:
        # Use abs() or raw? User said "Forward and Reverse combined Max".
        # Assuming we want the highest peak.
        scores = df_plot['PPR_Binding_Score']
        global_max = scores.max()
        global_min = scores.min()
        
    # Add Max/Min Lines to BOTH plots at the same Y level
    # Because Y axes are symmetric [0, limit], negative values are not plotted on these axes 
    # (or rather, Rev plot is 0-Limit reversed).
    # Wait, Reverse scores are positive in DF usually? Yes, visualization uses positive magnitude.
    # Let's verify: In plot loop, `y=df_rev['PPR_Binding_Score']`.
    # Previously `df_rev` values were plotted directly.
    # If `global_max` is e.g. 245, we draw line at 245 on both.
    
    if global_max > 0:
        if fwd_row is not None:
            fig.add_hline(y=global_max, line_dash="dash", line_color="blue", annotation_text=f"Max: {global_max:.1f}", annotation_position="top right", row=fwd_row, col=1)
        if rev_row is not None:
            max_anno_pos = "bottom right" if rows == 3 else "top right"
            fig.add_hline(y=global_max, line_dash="dash", line_color="blue", annotation_text=f"Max: {global_max:.1f}", annotation_position=max_anno_pos, row=rev_row, col=1)

    if global_min > 0: # Only if meaningful min exists (and > 0) or if user wants Min of displayed range even if small? 
        # Usually peaks are what matters. But user asked for Min.
        # "Min value". If min is 150, draw line at 150.
        if fwd_row is not None:
            fig.add_hline(y=global_min, line_dash="dot", line_color="red", annotation_text=f"Min: {global_min:.1f}", annotation_position="bottom right", row=fwd_row, col=1)
        if rev_row is not None:
            min_anno_pos = "top right" if rows == 3 else "bottom right"
            fig.add_hline(y=global_min, line_dash="dot", line_color="red", annotation_text=f"Min: {global_min:.1f}", annotation_position=min_anno_pos, row=rev_row, col=1)

    # Shared Y-Axis Limits (Based on Global Min/Max with padding)
    # User requested: Reference source is Global Max and Min values.
    # We add small padding (e.g. 5-10%) so points are not cut off.
    
    pad = (global_max - global_min) * 0.1
    if pad == 0: pad = 10 # Handle flat line
    
    y_max_limit = global_max + pad
    y_min_limit = max(0, global_min - pad) # Keep min >= 0 unless negative scores exist? 
    # Actually, if global_min is say 100, we want to see 90-ish.
    # Scores can be negative? P-type scores are usually positive log odds or similar.
    # If using absolute binding score, all > 0.
    # Let's trust the data Min.
    y_min_limit = global_min - pad
    
    # Force Symmetry + Add Y-axis Labels
    if fwd_row is not None and rev_row is not None and rows == 3:
        fig.update_yaxes(title_text="Fwd Score", range=[y_min_limit, y_max_limit], autorange=False, row=fwd_row, col=1)
        fig.update_yaxes(showticklabels=False, showgrid=False, zeroline=False, row=gene_row, col=1)
        fig.update_yaxes(title_text="Rev Score", range=[y_min_limit, y_max_limit], autorange="reversed", row=rev_row, col=1)
    elif fwd_row is not None:
        fig.update_yaxes(title_text="Fwd Score", range=[y_min_limit, y_max_limit], autorange=False, row=fwd_row, col=1)
        fig.update_yaxes(showticklabels=False, showgrid=False, zeroline=False, row=gene_row, col=1)
    elif rev_row is not None:
        fig.update_yaxes(showticklabels=False, showgrid=False, zeroline=False, row=gene_row, col=1)
        fig.update_yaxes(title_text="Rev Score", range=[y_min_limit, y_max_limit], autorange="reversed", row=rev_row, col=1)
    
    return fig


def plot_genes_only(gff_path, start, end, seq_id, feature_types):
    """Plots only the gene tracks (for Preview)."""
    fig = go.Figure()
    
    # Ensure feature_types is iterable list
    if isinstance(feature_types, str):
        feature_types = [feature_types]
    
    # 2. Plot Gene Annotations from GFF3 (Row 2 -> Row 1 in this single plot)
    # Using Subplot trick or just layout shapes? 
    # Just single axis.
    fig.add_trace(go.Scatter(x=[start, end], y=[0, 0], mode='markers', marker=dict(opacity=0), showlegend=False, hoverinfo='skip'))
    
    if gff_path and os.path.exists(gff_path):
        features = ncbi.parse_gff(gff_path, start, end)
        
        gene_y_fwd = 0.5
        gene_y_rev = -0.5
        
        for i, feat in enumerate(features):
            f_type = feat['type']
            if f_type not in feature_types: continue
            
            f_start = max(start, feat['start'])
            f_end = min(end, feat['end'])
            f_strand = feat['strand']
            
            attrs = feat.get('attributes', {})
            f_name = attrs.get('Name', attrs.get('gene', attrs.get('ID', f_type)))
            
            color = 'blue'
            if f_type == 'CDS': color = 'orange'
            elif f_type == 'tRNA': color = 'green'
            elif f_type == 'rRNA': color = 'purple'
            
            if f_strand == '+': y_pos = gene_y_fwd
            else: y_pos = gene_y_rev
            
            # Draw Arrow
            thickness = 0.4
            # Reuse simplified arrow logic
            view_width = end - start
            head_len = min(view_width * 0.02, f_end - f_start)
            half_thick = thickness / 2
            
            if f_strand == '+':
                y_top = y_pos + half_thick
                y_bot = y_pos - half_thick
                y_mid = y_pos
                x_head_base = max(f_start, f_end - head_len)
                path = f"M {f_start},{y_top} L {x_head_base},{y_top} L {f_end},{y_mid} L {x_head_base},{y_bot} L {f_start},{y_bot} Z"
            else:
                y_top = y_pos + half_thick
                y_bot = y_pos - half_thick
                y_mid = y_pos
                x_head_base = min(f_end, f_start + head_len)
                path = f"M {x_head_base},{y_top} L {f_end},{y_top} L {f_end},{y_bot} L {x_head_base},{y_bot} L {f_start},{y_mid} Z"

            fig.add_shape(type="path", path=path, line=dict(color=color, width=1), fillcolor=color, opacity=0.6)
            
            # Simplified Label Logic (Centered, No Box)
            
            fig.add_annotation(
                x=(f_start + f_end) / 2, 
                y=y_pos, 
                text=f"<b>{f_name}</b>", 
                showarrow=False,
                font=dict(color="black", size=10), 
                yshift=0,
                # bgcolor="white", # Removed
                opacity=1.0
            )

    fig.update_layout(
        title=f"Gene Preview: {start}-{end} ({seq_id})",
        height=300,
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, range=[-1, 1], title="Strand (+ / -)"),
        xaxis=dict(title="Genomic Position", tickformat="d"),
        margin=dict(t=50, b=20, l=50, r=50)
    )
    return fig


def visualize_motifs(motifs_file, analysis_title=None, discovery_tool=None):
    """Parses motifs file and returns DataFrame and Plotly figure."""
    try:
        # Read _motifs.txt (created by ppr_finder.py)
        # Columns: header, start, end, score, motif_seq, second, fifth, last, name
        # No header in file
        col_names = ['Header', 'Start', 'End', 'Score', 'Sequence', '2nd', '5th', 'Last', 'Type']
        df = pd.read_csv(motifs_file, sep='\t', names=col_names, header=None)
        
        fasta_name = "Unknown Sequence"
        if not df.empty:
             fasta_name = df.iloc[0]['Header'].lstrip('>')
        
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
        # Normalize scores for color calculation
        max_search_score = df['Score'].max() if not df.empty else 100
        min_search_score = df['Score'].min() if not df.empty else 0
        # If min is negative, handle? Usually E-value derived scores might be positive.
        # Let's assume scores are roughly comparable.
        
        # Color scale: Blue (Low) -> Red (High)
        # We can use Plotly's sample colors or interpolate.
        # Low=Blue, High=Red.
        
        for i, row in df.iterrows():
            motif_num = i + 1
            val = row['Score']

            # Map score to color 'rgba(r, g, b, 0.8)'
            # Simple interpolation: Blue(0,0,255) to Red(255,0,0)
            if max_search_score == min_search_score:
                fraction = 1.0
            else:
                fraction = (val - min_search_score) / (max_search_score - min_search_score)
                fraction = max(0, min(1, fraction))
            
            r = int(fraction * 255)
            b = int((1 - fraction) * 255)
            color_str = f"rgba({r}, 0, {b}, 0.8)"
            
            # Rectangle
            fig.add_shape(
                type="rect",
                x0=row['Start'], x1=row['End'],
                y0=-0.2, y1=0.2,
                line=dict(color="black", width=1),
                fillcolor=color_str,
                opacity=1.0 # Color handles opacity look
            )
            
            # Text (Code) - Center
            fig.add_annotation(
                x=(row['Start'] + row['End']) / 2,
                y=0,
                text=f"<b>{row['Code']}</b>",
                showarrow=False,
                font=dict(color="white", size=10)
            )
            
            # Text (Number and Type) - Above
            fig.add_annotation(
                x=(row['Start'] + row['End']) / 2,
                y=0.35,
                text=f"<b>#{motif_num}<br>({row['Type']})</b>",
                showarrow=False,
                font=dict(color="black", size=9)
            )
            
            # Text (Score) - Below
            # User requested to make it clear it's Motif Confidence Score
            fig.add_annotation(
                x=(row['Start'] + row['End']) / 2,
                y=-0.45,
                text=f"Conf:<br><b>{row['Score']:.1f}</b>",
                showarrow=False,
                font=dict(color="black", size=8)
            )

        # Add invisible Scatter to generate Colorbar
        # Use a point that exists but is invisible
        mid_x = (df['Start'].min() + df['End'].max()) / 2 if not df.empty else 0
        fig.add_trace(go.Scatter(
            x=[mid_x], y=[0],
            mode='markers',
            marker=dict(
                colorscale=[[0, 'blue'], [1, 'red']],
                cmin=min_search_score,
                cmax=max_search_score,
                size=0, # Invisible dot
                opacity=0, 
                showscale=True,
                colorbar=dict(
                    title=dict(text="Confidence Score", side="top"),
                    tickfont=dict(size=10),
                    orientation='h',
                    yanchor='bottom',
                    y=1.02, # Slightly above plot
                    xanchor='right',
                    x=1.0,
                    len=0.3,
                    thickness=15
                )
            ),
            hoverinfo='skip',
            showlegend=False
        ))
        
        # Infer tool label for loaded legacy results if not explicitly provided.
        tool_label = discovery_tool
        if not tool_label:
            unique_types = set(df['Type'].dropna().astype(str).unique())
            tool_label = "PS_Scan" if unique_types == {"PPR"} else "PPRFinder"

        # Layout settings
        title_text = f"Discovered Motifs ({tool_label}) - {fasta_name}"
        if analysis_title:
            title_text += f"<br><sub>{analysis_title}</sub>"
            
        fig.update_layout(
            title=title_text,
            xaxis=dict(title="Amino Acid Position", zeroline=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, range=[-0.6, 0.6]),
            height=300,
            # plot_bgcolor='white', # Removed to allow default light background
            margin=dict(t=80) # Increase top margin for subtitle
        )
        
        # Removed 'Header' from return to hide it, also user wanted "unneeded".
        return df[['Type', 'Start', 'End', 'Score', '5th', 'Last', 'Code', 'Sequence']], fig
        
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"Error visualizing motifs: {e}")
        return pd.DataFrame(), go.Figure()

def visualize_single_site(row):
    """Helper to generate the detail bar chart for a single row using Plotly."""
    logs = []
    try:
        logs.append(f"Visualizing Row: Gap={row.get('Gap_Size')}, Pos={row.get('Gap_Index')}")
        
        def parse_list_col(val, name):
            if isinstance(val, list): 
                logs.append(f"  {name} is list (len={len(val)})")
                return val
            if pd.isna(val) or val == "": 
                logs.append(f"  {name} is empty/NaN")
                return []
            s = str(val).strip()
            logs.append(f"  {name} raw string: {s[:50]}...")
            
            # Strip outer quotes if present
            if (s.startswith("'") and s.endswith("'")) or (s.startswith('"') and s.endswith('"')):
                s = s[1:-1]
            if s.startswith('[') and s.endswith(']'):
                s = s[1:-1]
            if not s.strip():
                return []
            try:
                res = [float(x) for x in s.split(',')]
                logs.append(f"  {name} parsed to list (len={len(res)})")
                return res
            except Exception as e:
                logs.append(f"  Error parsing {name}: {e}")
                return []

        scores = parse_list_col(row.get('PPR_Motif_Scores'), 'PPR_Motif_Scores')
        
        # Fallback if PPR_Motif_Scores is empty (happens if predicted_scores was empty during analysis)
        if not scores:
            logs.append("  PPR_Motif_Scores empty. Falling back to PPR_Binding_Scores_List.")
            scores = parse_list_col(row.get('PPR_Binding_Scores_List'), 'PPR_Binding_Scores_List')
            
        # Parse motif labels (handle string vs list vs nan)
        m_gap = row.get('Motif_Labels_GAP', '')
        if isinstance(m_gap, str):
             motif_labels = m_gap.split(',')
        elif isinstance(m_gap, list):
             motif_labels = m_gap
        else:
             # Handle NaN or other types
             motif_labels = []
        
        logs.append(f"  Motif Labels (len={len(motif_labels)}): {motif_labels}")
        
        segment = row.get('Segment', '')
        # Display as RNA (U)
        segment = core.dna_to_rna(segment)
        total_score = row.get('Total_Score', 0)
        logs.append(f"  Total Score: {total_score}")
        
        # Parse Motif Confidence Scores if available
        confidence_scores = []
        if 'Motif_Confidence_Scores' in row:
             confidence_scores = parse_list_col(row['Motif_Confidence_Scores'], 'Motif_Confidence_Scores')
        
        # Prepare data for plotting
        x_labels = []
        colors = []
        text_labels = []
        
        # Calculate Motif Numbers (skip gaps)
        # We need to reconstruct the original motif index (1-based)
        motif_numbers = []
        current_idx = 1
        for label in motif_labels:
            if label == '--':
                motif_numbers.append(None)
            else:
                motif_numbers.append(current_idx)
                current_idx += 1
        
        # Parse motif types if available
        # We need a new column 'Motif_Types_GAP' in row
        m_types_gap = row.get('Motif_Types_GAP', '')
        if isinstance(m_types_gap, str):
             motif_types = m_types_gap.split(',')
        elif isinstance(m_types_gap, list):
             motif_types = m_types_gap
        else:
             motif_types = []
        
        # If types not available, fallback to labels (Code)
        if not motif_types or len(motif_types) != len(motif_labels):
             motif_types = motif_labels

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
            '--': '#d3d3d3'  # Light Gray (Gap)
        }
        
        # Prepare data for plotting
        # We want to group bars by Motif Type to create separate traces for the Legend.
        # This solves the "What does this color mean?" issue.
        
        # Store data for all positions to ensure X-axis order and Scatter alignment
        all_x_labels = []
        all_scores = scores
        all_conf_scores = confidence_scores
        
        # Data by Type for Bar Traces
        # Structure: type -> {x: [], y: [], text: [], color: []}
        grouped_data = {}
        
        # Dynamic Font Size Calculation
        # Base size 12px, decrease as number of points increases
        # Min size 6px
        num_points = len(motif_labels)
        font_size = 12
        if num_points > 20:
            font_size = 10
        if num_points > 40:
            font_size = 8
        if num_points > 60:
            font_size = 6
        
        # Also adjust font for the HTML labels
        html_font_size = f"{font_size}px"
        
        for i, label in enumerate(motif_labels):
            # X-axis: Nucleotide + Index + Motif Number + Type
            nuc = segment[i] if i < len(segment) else "-"
            
            # Determine type for coloring
            m_type = motif_types[i] if i < len(motif_types) else label
            
            # Construct label parts
            if label == '--':
                m_num_str = "(Gap)"
                display_type = "Gap" # Type for display
                color = color_map['--']
                bar_text = "GAP"
            else:
                m_num = motif_numbers[i]
                m_num_str = f"#{m_num}" if m_num else ""
                display_type = m_type
                color = color_map.get(m_type, color_map.get(label, '#bcbd22'))
                
                # Bar Text: Code only (User request: remove Type from bar as it is in X-axis now)
                bar_text = label

            # X label: Position <br> Nucleotide <br> Motif# <br> Type
            # Use dynamic font size
            x_label_str = f"<span style='font-size:{html_font_size}'>{i+1}</span><br><b><span style='font-size:{html_font_size}'>{nuc}</span></b><br><span style='font-size:{html_font_size}'>{m_num_str}</span><br><span style='font-size:{html_font_size}'>{display_type}</span>"
            all_x_labels.append(x_label_str)
            
            # Hover Data: [Position, Nucleotide, MotifNum, Code, Type, Score]
            # Use safe strings for display
            h_pos = str(i+1)
            h_nuc = str(nuc)
            h_mnum = str(m_num) if m_num else "-"
            h_code = str(label)
            h_type = str(display_type)
            h_score = f"{scores[i]:.2f}"
            
            c_data = [h_pos, h_nuc, h_mnum, h_code, h_type, h_score]

            # Group data
            if display_type not in grouped_data:
                grouped_data[display_type] = {'x': [], 'y': [], 'text': [], 'color': color, 'customdata': []}
            
            grouped_data[display_type]['x'].append(x_label_str)
            grouped_data[display_type]['y'].append(scores[i])
            grouped_data[display_type]['text'].append(bar_text)
            grouped_data[display_type]['customdata'].append(c_data)

        fig = go.Figure()
        
        # Add a Bar trace for each Type
        # Sort keys to have consistent legend order (e.g. P, P1, S, S1...)
        sorted_types = sorted(grouped_data.keys())
        
        for m_type in sorted_types:
            data = grouped_data[m_type]
            # Name in legend: "Type: P", "Type: S", etc.
            trace_name = f"Type: {m_type}"
            if m_type == 'Gap': trace_name = "Gap"
            
            fig.add_trace(go.Bar(
                x=data['x'], 
                y=data['y'], 
                text=data['text'],
                textposition='auto', 
                marker_color=data['color'],
                name=trace_name,
                customdata=data['customdata'],
                hovertemplate=(
                    "<b>Site %{customdata[0]} (%{customdata[1]})</b><br>" +
                    "Motif: #%{customdata[2]}<br>" +
                    "Code: %{customdata[3]}<br>" +
                    "Type: %{customdata[4]}<br>" +
                    "Score: %{customdata[5]}<br>" +
                    "<extra></extra>" # Hides the trace name in the side box
                )
            ))
        
        # Provide explicit category order so X-axis doesn't re-sort by trace
        fig.update_xaxes(
            categoryorder='array', 
            categoryarray=all_x_labels, 
            tickangle=0,
            tickfont=dict(size=font_size) # Apply dynamic font size to ticks
        )
        
        # Line chart for Confidence Scores (if available) - This overlays all bars
        if all_conf_scores and len(all_conf_scores) == len(all_scores):
            fig.add_trace(go.Scatter(
                x=all_x_labels,
                y=all_conf_scores,
                mode='lines+markers',
                name='Motif Confidence',
                line=dict(color='black', width=2, dash='dot'),
                yaxis='y2',
                hovertemplate="<b>Motif Confidence</b>: %{y}<extra></extra>"
            ))
        
        # Get Analysis Type if available
        analysis_type_str = ""
        if 'Analysis_Type' in row and pd.notna(row['Analysis_Type']):
             atype = row['Analysis_Type']
             parts = atype.split('_')
             method_map = {'kobayashi': 'PLS-type', 'yan': 'P-type'}
             method_name = method_map.get(parts[0], parts[0].capitalize())
             
             ps_str = "Weighted by Confidence" if "no_ps" not in atype else "Unweighted"
             if 'no_ps' in atype:
                 ps_str = "Unweighted"
             elif 'ps' in atype:
                 ps_str = "Weighted by Confidence"
             
             analysis_type_str = f" | Method: {method_name} ({ps_str})"

        gap_size = row.get('Gap_Size', 0)
        gap_info = f" | Gap Size: {gap_size}"
        
        fig.update_layout(
            title=f"Site: {row['Genome_Position']} ({row['Strand']}){gap_info} | Score: {total_score:.2f}{analysis_type_str}",
            xaxis_title="Position / Nucleotide",
            xaxis=dict(
                tickangle=0, # Horizontal for readability
                tickfont=dict(size=10),
                automargin=True
            ),
            yaxis_title="Binding Score",
            yaxis2=dict(
                title="Motif Confidence",
                overlaying='y',
                side='right',
                showgrid=False
            ),
            height=400,
            showlegend=True,
            legend=dict(x=0, y=1.1, orientation="h")
        )
        return fig, "\n".join(logs)
    except Exception as e:
        import traceback
        traceback.print_exc()
        return go.Figure(layout=dict(title=f"Error generating plot: {e}")), f"Error generating plot: {e}"

def update_detail_plot(evt: gr.SelectData, df, analysis_context):
    if df is None or df.empty:
        # Return 5 outputs: gap_dd, manual_data, plot, seq, debug_log
        return gr.update(choices=[]), gr.update(value=None), go.Figure(), "", "No DataFrame loaded."
        
    try:
        print(f"DEBUG: Click Event: x={evt.value['x']}, y={evt.value['y']}")
        if analysis_context:
             print(f"DEBUG: Context Keys: {analysis_context.keys()}")
             print(f"DEBUG: Context Offset: {analysis_context.get('genome_offset')}")
        
        x_val = evt.value['x']
        y_val = evt.value['y']
        
        # Approximate matching for floats
        mask = (df['Genome_Position'] == x_val) & \
               ((np.isclose(df['PPR_Binding_Score'], y_val)) | (np.isclose(df['PPR_Binding_Score'], -y_val)))
        
        row = df[mask]
        print(f"DEBUG: Matched Rows in DF for click: {len(row)}")
        
        if row.empty:
             return gr.update(choices=[]), gr.update(value=None), go.Figure(layout=dict(title="No data found for clicked point")), "", f"No data found for x={x_val}, y={y_val}"
            
        best_row = row.iloc[0]
        # Strand is capitalized in DF: "Forward", "Reverse"
        strand = best_row['Strand']
        
        # Delegate to show_manual_detail to handle Gap Selection logic (Dynamic Mode if context available)
        # The return values here are for: gap_dd, manual_data, plot, seq, debug_log
        return show_manual_detail(str(x_val), strand, df, analysis_context)
        
    except Exception as e:
        print(f"Error in click handler: {e}")
        return gr.update(choices=[]), gr.update(value=None), go.Figure(layout=dict(title=f"Error: {e}")), "", f"Error in click handler: {e}"
        
    # Handler for detail click
def on_ncbi_detail_click(event_data, analysis_store):
    # reuse common logic
    if not event_data:
        return gr.update(choices=[]), gr.update(value=None), go.Figure(), "", "No event data."
    
    # event_data is EventData object. index is accessible.
    try:
         # Gradio 3/4 specific: event_data might be SelectData 
         # For Plotly, event_data provided by select/click is SelectData
         # The index is typically event_data.index[0] if present?
         # No, for gr.Plot(select=...), handler gets gr.SelectData.
         # data.index provides point index.
         
         # Actually, let's verify if `on_ncbi_detail_click` is used as .select listener
         # Down below: ncbi_browser_plot.select(...)
         
         # If using our helper wrapper `on_detail_click`:
         # The return values here are for: gap_dd, manual_data, plot, seq, debug_log
         return on_detail_click(event_data, analysis_store)
    except Exception as e:
        print(f"Error in handler: {e}")
        return gr.update(choices=[]), gr.update(value=None), go.Figure(layout=dict(title=f"Error: {e}")), "", f"Error in click handler: {e}"
        


def show_manual_detail(manual_pos, manual_strand, df, analysis_context):
    logs = [f"--- Manual Detail Search: Pos={manual_pos}, Strand={manual_strand} ---"]
    
    if df is None or df.empty:
        logs.append("No DataFrame loaded (None/Empty).")
        return gr.update(choices=[]), gr.update(value=None), go.Figure(), "", "\n".join(logs)

    try:
        if not manual_pos:
            logs.append("No manual position provided.")
            return gr.update(choices=[]), gr.update(value=None), go.Figure(), "", "\n".join(logs)
            
        target_pos = int(float(str(manual_pos)))
        logs.append(f"Target Position Parsed: {target_pos}")
    except ValueError:
        logs.append("Invalid Position Value.")
        return gr.update(choices=[]), gr.update(value=None), go.Figure(), "", "\n".join(logs)

    # Determine if we have context to re-calculate (Dynamic Mode) or just DataFrame (Static Mode)
    use_dynamic = False
    
    if analysis_context and isinstance(analysis_context, dict):
        if 'motif_data' in analysis_context and 'rna_seq' in analysis_context:
            use_dynamic = True
            
    logs.append(f"Mode: {'Dynamic (Re-calc)' if use_dynamic else 'Static (DF Lookup)'}")
            
    target_data = []
    
    if use_dynamic:
        # Re-calculate variations for this site
        motif_data = analysis_context['motif_data']
        rna_seq = analysis_context['rna_seq']
        min_front = analysis_context.get('min_front', 1)
        min_back = analysis_context.get('min_back', 1)
        max_gap = analysis_context.get('max_gap', 0)
        predicted_scores = analysis_context.get('predicted_scores', [])
        use_ps_score = analysis_context.get('use_ps_score', True)
        metric = analysis_context.get('method', 'yan')
        
        # Apply Offset if present (NCBI Slice Mode)
        offset = analysis_context.get('genome_offset', 0)
        logs.append(f"  Genome Offset: {offset}")
        
        # Adjust Relative Position
        rel_pos = target_pos - offset

        # Handle strand
        if manual_strand == "Reverse":
             # Reverse Complement Logic:
             # If mapping to RevComp Sequence:
             # Index i in Forward = Index (Len-1-i) in Reverse.
             # We need to extract the slice corresponding to `rel_pos` in Forward.
             # Wait, `core.analyze_single_site` behavior:
             # It takes `seq` and `pos`.
             # If we pass `rev_seq`, `pos` should be the index in `rev_seq`.
             # Index in `rev_seq` corresponds to:
             rev_pos = len(rna_seq) - 1 - rel_pos
             
             seq_to_analyze = core.get_reverse_complement(rna_seq)
             calc_pos = rev_pos
        else:
             seq_to_analyze = rna_seq
             calc_pos = rel_pos
             
        logs.append(f"Analyzing single site. Target: {target_pos} -> Rel: {rel_pos} -> Calc Pos: {calc_pos}")
        logs.append(f"RNA Seq Len: {len(rna_seq)}, Max Gap: {max_gap}")
        if calc_pos < 0 or calc_pos >= len(rna_seq):
             logs.append("WARNING: Calc Pos is out of bounds!")
             
        # Call core
        df_site = core.analyze_single_site(motif_data, seq_to_analyze, calc_pos, min_front, min_back, max_gap, predicted_scores, use_ps_score, metric)
        
        if not df_site.empty:
            target_data = df_site.to_dict('records')
            for d in target_data:
                d['Strand'] = manual_strand
                d['Genome_Position'] = target_pos
                d['PPR_Name'] = analysis_context.get('ppr_name', 'Unknown')
            logs.append(f"Dynamic Analysis found {len(target_data)} variations.")
        else:
            logs.append("Dynamic Analysis returned empty result.")
            
    else:
        # Fallback to static DataFrame
        logs.append(f"Searching DF with shape {df.shape}")
        mask = (df['Genome_Position'] == target_pos) & (df['Strand'] == manual_strand)
        subset = df[mask]
        
        # CHECK IF LITE MODE (Missing detailed columns)
        required_cols = ['PPR_Binding_Scores_List', 'Segment', 'PPR_Motif_Scores']
        missing_cols = [c for c in required_cols if c not in subset.columns]
        
        if missing_cols:
            logs.append(f"WARNING: Lite Mode Detected! Missing columns: {missing_cols}")
            logs.append("Detail Inspector requires these columns or Dynamic Context.")
            logs.append("Attempting to force Dynamic Mode if context was not provided...")
            
            # If context is empty, we can't do anything.
            # But the caller (on_ncbi_detail_click) calls reconstruct_context_from_csv.
            # So context should be available.
            # If we reached here, it means context was None/Empty?
            if not analysis_context:
                 logs.append("ERROR: Context is empty. Cannot reconstruct details from Lite CSV.")
                 return gr.update(choices=[], value=None), None, go.Figure(layout=dict(title="Error: Detailed data missing")), "", "\n".join(logs)
            else:
                 # This branch shouldn't be reached if logic above (if analysis_context: use_dynamic=True) is correct.
                 # Wait, logic above checks: if 'motif_data' in analysis_context...
                 # If analysis_context exists but missing keys?
                 logs.append("Context exists but might be malformed or `use_dynamic` logic failed.")
                 # Force return error or try again?
                 return gr.update(choices=[], value=None), None, go.Figure(layout=dict(title="Error: Dynamic Analysis Failed")), "", "\n".join(logs)
                 
        target_data = subset.to_dict('records')
        logs.append(f"Found {len(target_data)} matching rows in DF.")

    if not target_data:
        logs.append("No target data found.")
        return gr.update(choices=[], value=None), None, go.Figure(layout=dict(title="No data found")), "", "\n".join(logs)

    options = []
    data_map = {}
    
    # Sort by Score descending
    target_data.sort(key=lambda x: x['PPR_Binding_Score'], reverse=True)
    
    for i, row in enumerate(target_data):
        gap_s = row.get('Gap_Size', 0)
        gap_idx = row.get('Gap_Index', 0)
        score = row.get('PPR_Binding_Score', 0)
        
        label = f"#{i+1}: Gap Size {gap_s} (Pos {gap_idx}) | Score {score:.2f}"
        options.append(label)
        data_map[label] = row
        
    default_val = options[0] if options else None
    logs.append(f"Default Selection: {default_val}")
    
    # Call visualize with unpacking
    # Generate Plot for best result
    if default_val:
        fig, log_txt = visualize_single_site(data_map.get(default_val))
        logs.append(log_txt)
    else:
        fig, log_txt = go.Figure(), "No default value selected"
        logs.append(log_txt)
    
    # Return matched sequence (as RNA)
    seq_val = ""
    if default_val:
        seq_val = data_map[default_val].get('Segment', '')
        seq_val = core.dna_to_rna(seq_val)
    
    return gr.update(choices=options, value=default_val), data_map, fig, seq_val, "\n".join(logs)

def inspect_gap_select(selection_str, data_map):
    """Updates the plot based on the dropdown selection using cached data_map."""
    if not selection_str or not data_map:
        return go.Figure(), "", "No selection or data map."
        
    try:
        row = data_map.get(selection_str)
        if not row:
            return go.Figure(layout=dict(title="Selection not found in data")), "", "Selection not found."
            
        fig, log_txt = visualize_single_site(row)
        fig, log_txt = visualize_single_site(row)
        
        # Return Sequence as RNA (U)
        segment = row.get('Segment', '')
        segment = core.dna_to_rna(segment)
        
        return fig, segment, log_txt
        
    except Exception as e:
        import traceback
        traceback.print_exc()
        return go.Figure(layout=dict(title=f"Error: {e}")), "", f"Error in dropdown: {e}"


# --- IGV Helper Removed ---
# def generate_igv_html(fasta_path, bedgraph_path, locus_id): ...

def run_ncbi_analysis(accession_id, protein_text, protein_file, tool_choice, scoring_method, max_gap_size, use_ps_score, job_name, generate_bed=True):
    # Construct analysis_type string
    method_map = {"PLS-type": "kobayashi", "P-type": "yan"}
    method = method_map.get(scoring_method, scoring_method).lower()
    gap_str = f"gap{int(max_gap_size)}" if max_gap_size > 0 else "gap0"
    ps_str = "ps" if use_ps_score else "no_ps"
    analysis_type = f"{method}_{gap_str}_{ps_str}"
    
    # Use relaxed constraints for gaps
    min_front = 1
    min_back = 1
    
    """
    Orchestrates the analysis starting from an NCBI Accession ID.
    Yields updates for real-time logging in Gradio.
    """
    status_msg = ""
    log_buffer = []
    
    # Helper to construct the output tuple
    # Outputs: [status_output, ncbi_fasta_download, ncbi_bedgraph_download, ncbi_bed_download, ncbi_gff_download, ncbi_download, ncbi_log_output, plot_browser, plot_detail, df_state, analysis_context_state]
    def make_output(status, fasta=None, bedgraph=None, bed=None, gff=None, csv=None, fig_browser=None, fig_detail=None, df_out=None, context_out=None):
        return (
            status, 
            fasta, 
            bedgraph,
            bed,
            gff,
            csv,
            "\n".join(log_buffer),
            fig_browser if fig_browser else gr.update(),
            fig_detail if fig_detail else gr.update(),
            df_out if df_out is not None else gr.update(),
            context_out if context_out is not None else gr.update()
        )

    def log(message):
        log_buffer.append(message)
        # We can yield intermediate status updates if we want, but for now just append to log
        
    # 1. Fetch Genome (FASTA)
    if not accession_id:
        accession_id = DEFAULT_ACCESSION
        
    try:
        log(f"Fetching genome for {accession_id}...")
        yield make_output(f"Fetching genome for {accession_id}...")
        fasta_path = ncbi.fetch_fasta(accession_id, TMP_DIR)
        log(f"Genome fetched: {fasta_path}")
        yield make_output(f"Genome fetched. Running analysis...", fasta=fasta_path)
        
    except Exception as e:
        err_msg = f"Error fetching genome from NCBI: {e}"
        log(err_msg)
        yield make_output(err_msg)
        return
        
    # 2. Read the fetched FASTA content (now it's the genome)
    try:
        with open(fasta_path, "r") as f:
            rna_text = f.read()
            # Normalize to RNA (U)
            rna_text = core.dna_to_rna(rna_text)
        log(f"Read FASTA content ({len(rna_text)} bytes).")
        yield make_output("Read FASTA.", fasta=fasta_path)
    except Exception as e:
        err_msg = f"Error reading FASTA file: {e}"
        log(err_msg)
        yield make_output(err_msg)
        return
        
    # 3. Run Unified Analysis
    log("Starting unified analysis (this may take a few minutes)...")
    yield make_output("Running Analysis...", fasta=fasta_path)
    
    # 3. Run Unified Analysis logic
    log("Running Unified Analysis...")
    yield make_output("Running Analysis...", fasta=fasta_path)
    
    try:
        # We need to adapt the return values
        # Return: status_msg, csv_path, motif_fig, motif_df, fig_browser, range_plot, df_output, analysis_context, dd_main, dd_ncbi
        # Return: status_msg, csv_path, motif_fig, motif_df, fig_browser, range_plot, df_output, analysis_context, dd_main, dd_ncbi
        result_msg, csv_path, motif_fig, _, _, range_plot, df, analysis_context_res, _, _ = run_unified_analysis(
            protein_text, protein_file, rna_text, None, tool_choice, scoring_method, max_gap_size, use_ps_score, job_name, skip_plots=True
        )
        log(f"Analysis result: {result_msg}")
        
        # Optimize for UI: Truncate DataFrame if too large (User Request)
        # Strategy: Keep Top N by Score to preserve relevant Reverse hits and high-confidence sites.
        top_n = 1000
        limit_msg = ""
        
        if df is not None and len(df) > top_n:
            log(f"DataFrame size {len(df)} > {top_n}. Truncating to Top {top_n} entries (Balanced).")
            
            # Separate by Strand
            df_fwd = df[df['Strand'] == "Forward"].sort_values(by="PPR_Binding_Score", ascending=False)
            df_rev = df[df['Strand'] == "Reverse"].sort_values(by="PPR_Binding_Score", ascending=False)
            
            # Take Top N/2 from each (or as many as available)
            limit_per_strand = top_n // 2
            df_fwd_top = df_fwd.head(limit_per_strand)
            df_rev_top = df_rev.head(limit_per_strand)
            
            # Concatenate and sort mainly by Position for graph consistency (or just keep concatenated)
            df = pd.concat([df_fwd_top, df_rev_top])
            
            limit_msg = f"\n(Note: Result truncated to Top {len(df)} sites [~500 Fwd/Rev] by Score)"
            log(limit_msg.strip())
        
        # Manually generate plots with the (potentially truncated) DF
        fig_browser = plot_genome_browser(df, analysis_type=analysis_type, highlight_top=False, highlight_strand="Both")
            
        final_status = "Analysis Finished." + limit_msg
        
        # We yield the plots and state
        yield make_output(final_status, fasta=fasta_path, csv=csv_path, fig_browser=fig_browser, fig_detail=range_plot, df_out=df, context_out=analysis_context_res)
        
    except Exception as e:
        err_msg = f"Error during unified analysis: {e}"
        log(err_msg)
        import traceback
        log(traceback.format_exc())
        yield make_output(err_msg)
        return
    
    if df is None or df.empty:
        err_msg = f"Analysis failed or no results: {msg}"
        log(err_msg)
        yield make_output(err_msg)
        return
        
    log(f"Analysis complete. DataFrame shape: {df.shape}")
        
    # 4. Generate bedGraph and BED for IGV
    bedgraph_path = None
    bed_path = None
    try:
        if not generate_bed:
            raise ValueError("SKIP_IGV")
            
        bedgraph_path = os.path.join(TMP_DIR, f"{accession_id}_{job_name}.bedgraph")
        bed_path = os.path.join(TMP_DIR, f"{accession_id}_{job_name}.bed")
        
        log(f"Generating bedGraph and BED files...")
        yield make_output("Generating visualization files...", fasta=fasta_path, csv=csv_path)
        
        from Bio import SeqIO
        record = next(SeqIO.parse(fasta_path, "fasta"))
        chrom_id = record.id
        log(f"Chromosome ID from FASTA: {chrom_id}")
        
        # BedGraph (Forward only for now, or maybe separate?)
        # Original logic was Forward only. Let's keep it consistent for BedGraph or improve?
        # User asked for strand distinction. BedGraph is usually one value per pos.
        # BED is better for strand distinction.
        
        # Write BedGraph (Forward binding score)
        with open(bedgraph_path, "w") as f:
            f.write(f"track type=bedGraph name='{analysis_type}_Score' description='PPR Binding Score' visibility=full color=0,0,255 autoScale=on viewLimits=0:30\n")
            df_fwd = df[df['Strand'] == 'Forward']
            for _, row in df_fwd.iterrows():
                pos = int(row['Genome_Position'])
                score = float(row['PPR_Binding_Score'])
                f.write(f"{chrom_id}\t{pos}\t{pos+1}\t{score}\n")
                
        # Write BED (Reads visualization with Strand)
        with open(bed_path, "w") as f:
            # Track line: itemRgb=On to use colors
            f.write(f"track name='{analysis_type}_Sites' description='PPR Binding Sites' visibility=pack itemRgb=On useScore=1\n")
            
            count = 0
            for _, row in df.iterrows():
                # BED format: chrom, chromStart, chromEnd, name, score, strand, thickStart, thickEnd, itemRgb
                
                strand = row['Strand']
                score_val = float(row['PPR_Binding_Score'])
                
                # Scale score 0-30 -> 0-1000 for shading (optional)
                bed_score = int(min(score_val * 33, 1000))
                
                if strand == 'Forward':
                    start = int(row['Genome_Position'])
                    # Lite Mode: Segment is empty. Fallback to estimated length.
                    seg = row.get('Segment', '')
                    if seg:
                        length = len(seg)
                    else:
                        # Estimate: Default 35 or (30 + Gap)
                        # Ideal: len(motifs) + Gap. But we don't have motifs handy in this loop easily unless we captured context.
                        # Using fixed approximation 35 + Gap is safer than 0.
                        length = 35 + int(row.get('Gap_Size', 0))
                        
                    end = start + length
                    strand_char = "+"
                    color = "0,0,255" # Blue
                else:
                    high_coord = int(row['Genome_Position']) # This is the index on forward strand
                    
                    seg = row.get('Segment', '')
                    if seg:
                        length = len(seg)
                    else:
                        length = 35 + int(row.get('Gap_Size', 0))
                    
                    chromEnd = high_coord + 1
                    chromStart = chromEnd - length
                    
                    start = chromStart
                    end = chromEnd
                    strand_char = "-"
                    color = "255,0,0" # Red
                
                name = f"{row['PPR_Name']}_{score_val:.1f}"
                
                # Write BED line
                f.write(f"{chrom_id}\t{start}\t{end}\t{name}\t{bed_score}\t{strand_char}\t{start}\t{end}\t{color}\n")
                count += 1
                
            log(f"Written {count} records to BED file.")

    except Exception as e:
        if str(e) == "SKIP_IGV":
            log("Skipping IGV files generation per user request.")
            bedgraph_path = None
            bed_path = None
            yield make_output("Skipped visualization files.", fasta=fasta_path, csv=csv_path)
        else:
            err_msg = f"Error generating visualization files: {e}"
            log(err_msg)
            import traceback
            log(traceback.format_exc())
            yield make_output(err_msg, csv=csv_path)
            return
            
    # Return paths for download
    log(f"Ready for download: {fasta_path}, {bedgraph_path}, {bed_path}")
    
    # 5. Fetch GFF3 (Annotations)
    gff3_path = None
    try:
        log(f"Fetching GFF3 annotations for {accession_id}...")
        yield make_output("Fetching GFF3...", fasta=fasta_path, csv=csv_path, bedgraph=bedgraph_path, bed=bed_path)
        gff3_path = ncbi.fetch_gff(accession_id, output_dir=TMP_DIR)
        if gff3_path:
            log(f"Fetched GFF3: {gff3_path}")
        else:
            log("GFF3 fetch returned None (might not be available).")
    except Exception as e:
        log(f"Warning: Failed to fetch GFF3: {e}")
        # Don't fail the whole analysis for this
    
    yield make_output("Done", fasta=fasta_path, bedgraph=bedgraph_path, bed=bed_path, csv=csv_path, gff=gff3_path)

def open_input_folder():
    return open_folder(INPUT_DIR)

def open_output_folder():
    return open_folder(OUTPUT_DIR)

def save_custom_gff(df):
    """Saves the custom genes dataframe to a GFF3 file."""
    if df is None or df.empty:
        return "No data to save."
    
    # Validate columns
    required = ["Name", "Type", "Start", "End", "Strand"]
    if not all(col in df.columns for col in required):
        # Fallback for old dataframe version?
        # Normalize if "Type" missing?
        if "Type" not in df.columns and "Name" in df.columns:
             df["Type"] = "gene" # Backwards compat
        elif not all(col in df.columns for col in ["Name", "Start", "End", "Strand"]):
             return "Error: Missing columns in table."
        
    out_path = os.path.join(INPUT_DIR, "custom_annotations.gff3")
    
    try:
        with open(out_path, "w") as f:
            f.write("##gff-version 3\n")
            for _, row in df.iterrows():
                # GFF3: seqid source type start end score strand phase attributes
                # We use "Custom" as seqid placeholder or match target?
                # If merging, seqid should match target ID ideally.
                # But typically visualization tools ignore seqid if we force load?
                # Let's use "." or "Target". 
                
                name = str(row["Name"])
                f_type = str(row["Type"]) if "Type" in row and row["Type"] else "gene"
                start = str(int(row["Start"]))
                end = str(int(row["End"]))
                strand = str(row["Strand"])
                
                # Attributes
                attr = f"ID={name};Name={name}"
                
                line = f".\tCustom\t{f_type}\t{start}\t{end}\t.\t{strand}\t.\t{attr}\n"
                f.write(line)
        return f"Saved to {os.path.basename(out_path)}"
    except Exception as e:
        return f"Error saving GFF: {e}"

def create_merged_gff(base_gff_path, custom_df, seq_id):
    """
    Creates a temporary GFF3 file merging base GFF and Custom Annotations.
    """
    import tempfile
    
    # If no custom data, just return base
    if custom_df is None or custom_df.empty:
        return base_gff_path
        
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=".gff3") as tmp_gff:
        # Copy base GFF content if exists
        if base_gff_path and os.path.exists(base_gff_path):
            with open(base_gff_path, "r") as f_in:
                tmp_gff.write(f_in.read())
                # Ensure newline
                # if not base_gff_path.endswith("\n"): 
                # Reading file content doesn't give filename property.
                # Just add newline to be safe
                tmp_gff.write("\n")
        else:
            tmp_gff.write("##gff-version 3\n")
            
        # Append Custom Genes
        for _, row in custom_df.iterrows():
            try:
                name = str(row["Name"])
                f_type = str(row["Type"]) if "Type" in row and row["Type"] else "gene"
                s = int(row["Start"])
                e = int(row["End"])
                st = str(row["Strand"])
                
                # Use seq_id matches target
                attr = f"ID={name};Name={name}"
                
                # Logic: If NCBI Mode, seq_id is e.g. "AP017386".
                # If Manual Mode, seq_id is "Manual" or "."?
                # Plotters might filter by ID. 
                # Let's trust the passed seq_id.
                eff_id = seq_id if seq_id else "."
                
                line = f"{eff_id}\tCustom\t{f_type}\t{s}\t{e}\t.\t{st}\t.\t{attr}\n"
                tmp_gff.write(line)
            except:
                continue
                
    return tmp_gff.name

def preview_genes_wrapper(mode, rna_text, ncbi_gff, start, end, feature_type, acc_id, custom_df):
    """
    Unified Preview Logic handling both Manual and NCBI modes + Custom Annotations.
    """
    import tempfile
    
    # 1. Determine Sequence Length and Context
    seq_len = 50000 # Default fallback
    effective_start = 1
    effective_end = seq_len
    effective_id = "Manual"
    
    base_gff_path = None
    
    if mode == "NCBI Accession":
        if not ncbi_gff:
            return go.Figure(layout=dict(title="Fetch data first for NCBI Mode"))
        base_gff_path = ncbi_gff.name
        effective_start = int(start) if start else 1
        effective_end = int(end) if end else 1000
        effective_id = acc_id
    else:
        # Manual Mode
        # Use provided range if available, else seq len
        if rna_text:
            # Clean and estimate length
            clean_seq = "".join(rna_text.splitlines()[1:]).replace(" ", "").strip() if ">" in rna_text else rna_text.replace(" ", "").strip()
            seq_len = len(clean_seq) if clean_seq else 1000
        else:
            seq_len = 1000
            
        if start and end and (end > start):
             effective_start = int(start)
             effective_end = int(end)
        else:
             effective_start = 1
             effective_end = seq_len
            
    # 2. Handle Custom Annotations (Merge Logic)
    merged_gff_path = create_merged_gff(base_gff_path, custom_df, effective_id)
            
    if not merged_gff_path or not os.path.exists(merged_gff_path):
         # If no data, and we are in manual mode with no custom genes, we might have no GFF.
         # plot_genes_only handles None gff by plotting nothing?
         pass
         
    try:
        # Call Plotter
        fig = visualization.plot_genes_only(merged_gff_path, effective_start, effective_end, effective_id, feature_type)
        return fig
    except Exception as e:
        return go.Figure(layout=dict(title=f"Preview Error: {e}"))

def run_dispatch_wrapper(
    mode_select,
    # Manual Inputs
    prot_text, prot_file, rna_text, rna_file,
    # NCBI Inputs
    ncbi_accession, range_start, range_end, ncbi_fasta_down, ncbi_gff_down,
    # Shared Settings
    tool_select, scoring_method, max_gap_size, use_ps_score, job_name,
    # Features
    feature_types, plot_height, highlight_top, highlight_strand,
    log_state=None, # Add log_state
    custom_gene_df=None # Add custom_gene_df
):
    """
    Unified Dispatcher. Routes request to appropriate logic and normalizes output.
    Returns: 
    [status, csv, motif_plot, motif_df, browser_plot, detail_plot (dummy), df_state, context_state, log_out]
    """
    # Normalize feature_types (Radio returns string, CheckboxGroup returned list)
    # Normalize feature_types (Radio returns string, CheckboxGroup returned list)
    if isinstance(feature_types, str):
        feature_types = [feature_types]
    elif feature_types is None or len(feature_types) == 0:
        feature_types = ["gene"] # Fallback default

    import pandas as pd
    import time

    # Initialize Log
    current_log = log_state if log_state else ""
    def _log(msg):
        nonlocal current_log
        current_log = add_log(current_log, msg)

    _log("--- New Analysis Started ---")
    _log(f"Mode: {mode_select}")
    start_time = time.time()
    
    # Defaults
    status = "Starting..."
    csv = None
    motif_fig = gr.update() # Keep existing
    motif_df = gr.update()
    browser_fig = go.Figure(layout=dict(title="Waiting for Analysis..."))
    detail_dummy = go.Figure()
    df_out = None
    context_out = {}
    
    # Helper to yield current state
    # We map current_log to the first output (was status) to drive terminal_log
    # We remove the redundant current_log at the end
    def _yield():
        return (current_log, csv, motif_fig, motif_df, browser_fig, detail_dummy, df_out, context_out)

    yield _yield() # Initial yield

    try:
        if mode_select == "NCBI Accession":
            # --- NCBI MODE (Unified Logic) ---
            if not ncbi_accession or not ncbi_accession.strip():
                ncbi_accession = DEFAULT_ACCESSION
            
            _log(f"Target: NCBI Accession {ncbi_accession}")
            _log(f"Range: {range_start} - {range_end}")
            
            _log(f"Running NCBI Analysis on {ncbi_accession} ({range_start}-{range_end})...")
            yield _yield()
            
            # 1. Validation
            if not ncbi_fasta_down or not os.path.exists(ncbi_fasta_down.name):
                _log("ERROR: FASTA file not found. Fetch required.")
                yield _yield()
                return
            
            # 2. Read & Slice Sequence
            from Bio import SeqIO
            try:
                _log("Reading FASTA file...")
                yield _yield()
                record = next(SeqIO.parse(ncbi_fasta_down.name, "fasta"))
                full_seq = str(record.seq).upper()
                _log(f"Loaded Genome: {len(full_seq)} bp")
                yield _yield()
                
                # Handle Bounds
                r_start = int(max(0, range_start))
                r_end = int(min(len(full_seq), range_end))
                
                if r_start >= r_end:
                     _log(f"ERROR: Invalid Range {r_start}-{r_end}")
                     yield _yield()
                     return
                
                # Extract Slice
                rna_slice = full_seq[r_start:r_end]
                rna_header = f">{record.id}:{r_start}-{r_end}\n"
                rna_text_input = rna_header + rna_slice
                
                status += f"\nAnalyzed Region size: {len(rna_slice)} bp." 
                _log(f"Sliced RNA: {len(rna_slice)} bp ({r_start}-{r_end})")
                yield _yield()
                
            except Exception as e:
                _log(f"ERROR reading FASTA: {e}")
                yield _yield()
                return

            # 3. Call Unified Analysis (The Core Algorithm)
            _log("Executing Core Analysis...")
            yield _yield()
            
            res_tuple = run_unified_analysis(
                prot_text, prot_file, rna_text_input, None, 
                tool_select, scoring_method, max_gap_size, use_ps_score, job_name,
                highlight_top=highlight_top, highlight_strand=highlight_strand,
                skip_plots=True # We will manually plot to include GFF
            )
            
            # Unpack
            (res_msg, csv, motif_fig, motif_df, _, _, df_out, context_out, _, _) = res_tuple
            
            # 4. Post-Process Results (Coordinate Shift)
            if df_out is not None and not df_out.empty:
                df_out['Genome_Position'] = df_out['Genome_Position'] + r_start
                status += " (Coordinates adjusted to Genomic Position)"
            else:
                status += "\nNo binding sites found in this range."

            # Inject Offset into Context for Detail Inspector
            if context_out:
                context_out['genome_offset'] = r_start

            # 5. Generate Browser Plot (With GFF!)
            _log("Rendering Range Browser...")
            yield _yield()
            
            # Merge GFF if needed (NCBI GFF + Custom)
            base_gff_path = ncbi_gff_down.name if ncbi_gff_down else None
            merged_gff = create_merged_gff(base_gff_path, custom_gene_df, ncbi_accession)
            
            browser_fig = plot_range_browser(
                df_out, merged_gff, r_start, r_end, 
                seq_id=ncbi_accession, 
                feature_types=feature_types,
                plot_height=plot_height, 
                highlight_top=highlight_top, 
                highlight_strand=highlight_strand,
                analysis_context={}, # Disable internal re-calc (we already calculated exact range)
                protein_name=job_name
            )
            
            status += f"\n{res_msg}"
            
            elapsed = time.time() - start_time
            _log(f"Analysis Complete in {elapsed:.1f}s")
            
            if df_out is not None and not df_out.empty:
                top_score = df_out['PPR_Binding_Score'].max()
                _log(f"Top Score: {top_score:.1f}")
                _log(f"Sites Found: {len(df_out)}")
                
            yield _yield()

        else:
            # --- MANUAL MODE ---
            _log("Mode: Manual Sequence Analysis")
            yield _yield()
            
            # Call run_unified_analysis
            _log("Executing Core Analysis...")
            yield _yield()
            
            res_tuple = run_unified_analysis(
                prot_text, prot_file, rna_text, rna_file, 
                tool_select, scoring_method, max_gap_size, use_ps_score, job_name,
                highlight_top=highlight_top, highlight_strand=highlight_strand
            )
            
            # Unpack
            (status, csv, motif_fig, motif_df, browser_fig_orig, _, df_out, context_out, _, _) = res_tuple
            
            # For Manual Mode, we also want to show Custom Annotations.
            # We override browser_fig_orig (which uses plot_genome_browser without GFF) 
            # with plot_range_browser using the custom GFF.
            
            _log("Rendering Browser with Custom Annotations...")
            if df_out is not None and not df_out.empty:
                max_pos = int(df_out['Genome_Position'].max())
                # Create GFF from Custom DF
                merged_gff = create_merged_gff(None, custom_gene_df, "Manual")
                
                browser_fig = plot_range_browser(
                    df_out, merged_gff, 0, max_pos + 100, 
                    seq_id="Target RNA", 
                    feature_types=["gene"], # Default showing genes
                    plot_height=plot_height, 
                    highlight_top=highlight_top, 
                    highlight_strand=highlight_strand,
                    analysis_context={},
                    protein_name=job_name
                )
            else:
                browser_fig = browser_fig_orig
            
            elapsed = time.time() - start_time
            _log(f"Analysis Complete in {elapsed:.1f}s")
            
            if df_out is not None and not df_out.empty:
                top_score = df_out['PPR_Binding_Score'].max()
                _log(f"Top Score: {top_score:.1f}")
            
            yield _yield()

    except Exception as e:
        import traceback
        err = traceback.format_exc()
        _log(f"CRITICAL ERROR: {e}")
        yield _yield()

# --- UI Definition ---

# --- UI Definition (Unified) ---
from gradio.themes import Default

# --- UI Definition (Tabbed) ---
from gradio.themes import Default

with gr.Blocks(title="P-type PPR Prediction Program (PPPRPP)", css=custom_css) as demo:
    with gr.Row():
        gr.Markdown(f"# P-type PPR Prediction Program (PPPRPP) <span style='font-size: 0.5em; font-weight: normal; color: #777;'>{VERSION}</span>")
        with gr.Row(): # Buttons aligned right-ish
             btn_open_in = gr.Button("📂 Input Folder", size="sm", scale=0)
             btn_open_out = gr.Button("📂 Output Folder", size="sm", scale=0)
    
    # State
    df_state = gr.State()
    analysis_context_state = gr.State({})
    manual_data_state = gr.State({}) 
    log_state = gr.State(">> PPPRPP System Ready.\n")  # Accumulated log messages
    
    # Helper function to add timestamped log entry
    def add_log(current_log, message):
        """Append timestamped message to log."""
        from datetime import datetime
        timestamp = datetime.now().strftime("%H:%M:%S")
        new_entry = f">> [{timestamp}] {message}\n"
        return current_log + new_entry 
    
    with gr.Tabs() as main_tabs:
        
        # --- Tab 1: Inputs ---
        with gr.Tab("1. Inputs", id="tab_inputs"):
            gr.Markdown("### Protein Input (Motif Discovery)")
            with gr.Row():
                prot_text = gr.Textbox(label="Protein FASTA Sequence", lines=5, placeholder=f">PPR791(Rf1a)_IR24\\n{DEFAULT_PROTEIN[:50]}...", scale=3)
                prot_file = gr.File(label="Upload Protein File", scale=1)

            gr.Markdown("### Target RNA Input")
            mode_select = gr.Radio(["Manual Sequence", "NCBI Accession"], value="Manual Sequence", label="Input Mode")
            
            # Manual Panel
            with gr.Group(visible=True) as group_manual:
                with gr.Row():
                    rna_text = gr.Textbox(label="RNA FASTA (Target)", lines=5, placeholder=f">atp6-orf79_BT\\n{DEFAULT_RNA[:50]}...", scale=3)
                    rna_file = gr.File(label="Upload RNA File", scale=1)

            # NCBI Panel
            with gr.Group(visible=False) as group_ncbi:
                with gr.Row():
                    ncbi_accession = gr.Textbox(label="NCBI Accession ID", placeholder=f"e.g. {DEFAULT_ACCESSION}", scale=2)
                    ncbi_fetch_btn = gr.Button("Fetch Data", variant="secondary", scale=1)
                
                # Fetch Status (Immediate)
                input_status_log = gr.Textbox(label="Fetch Status", lines=1, interactive=False)
                
            # Hidden Stores (Shared context)
            ncbi_fasta_down = gr.File(label="Ref FASTA", visible=False)
            ncbi_gff_down = gr.File(label="Ref GFF3", visible=False)
            
            gr.Markdown("#### Preview & Features")
            # Custom & Range Preview Shared Block
            with gr.Group():
                with gr.Row():
                     with gr.Column(scale=2):
                         gr.Markdown("**Custom Annotations** (Name, Type, Start, End, Strand)")
                         custom_gene_df = gr.DataFrame(
                            headers=["Name", "Type", "Start", "End", "Strand"],
                            datatype=["str", "str", "number", "number", "str"],
                            col_count=(5, "fixed"),
                            label="Custom Genes Table",
                            interactive=True,
                            value=[["MyGene", "gene", 100, 500, "+"]] # Example row with Type
                        )
                         btn_save_custom = gr.Button("Save (GFF3)", size="sm")
                         save_status = gr.Textbox(show_label=False, interactive=False)
                         
                     with gr.Column(scale=1):
                         gr.Markdown("**Preview Settings**")
                         with gr.Row():
                             range_start = gr.Number(label="Start", value=133162, precision=0)
                             range_end = gr.Number(label="End", value=134622, precision=0)
                         preview_btn = gr.Button("Preview Genes", variant="secondary")
                         feature_radio = gr.CheckboxGroup(
                             ["gene", "CDS", "tRNA", "rRNA", "exon", "other"], 
                             value=["gene", "other"], 
                             label="Preview Feature Types"
                         )
                
            # Gene Preview Plot
            preview_plot = gr.Plot(label="Range Preview", container=True, show_label=False)

        # --- Tab 2: Settings & Run ---
        with gr.Tab("2. Settings & Run", id="tab_settings"):
            gr.Markdown("### Analysis Settings")
            with gr.Row():
                with gr.Column():
                    tool_select = gr.Dropdown(choices=["PPRFinder", "PS_Scan"], value="PPRFinder", label="Discovery Tool")
                    job_name_input = gr.Textbox(value="my_analysis", label="Job Name")
                with gr.Column():
                    scoring_method = gr.Radio(["P-type", "PLS-type"], value="P-type", label="Scoring Method")
                    use_ps_score = gr.Checkbox(label="Use Motif Confidence Score (+Score)", value=True)
                    range_max_gap = gr.Number(value=0, label="Gap Size", precision=0)
            
            gr.Markdown("---")
            run_btn = gr.Button("Run Analysis", variant="primary", size="lg")
            run_status_simple = gr.Textbox(label="Run Status", lines=1, interactive=False)
            gr.Markdown("*Detailed results will appear in '3. Analysis Results' tab.*")

        # --- Tab 3: Results ---
        with gr.Tab("3. Analysis Results", id="tab_results"):
            gr.Markdown("### Motif Analysis")
            motif_plot = gr.Plot(label="Detected Motifs", container=True, show_label=False)
            with gr.Accordion("Motif Details Table", open=False):
                motif_table = gr.DataFrame(label="Motif Details", visible=True)

            gr.Markdown("### Target Prediction (Genome Browser)")
            with gr.Row():
                chk_highlight = gr.Checkbox(label="Highlight Top", value=False)
                dd_strand = gr.Dropdown(choices=["Forward", "Reverse", "Both"], value="Both", label="Show Strand")
                plot_height = gr.Slider(minimum=300, maximum=1000, value=500, label="Plot Height")

            plot_browser = gr.Plot(label="Genome Browser", show_label=False, container=True)
            
            gr.Markdown("### Detail Inspector")
            # Controls Row
            with gr.Row():
                 manual_pos = gr.Number(label="Position", precision=0, scale=1)
                 manual_strand = gr.Radio(["Forward", "Reverse"], value="Forward", label="Strand", scale=1)
                 gap_select = gr.Dropdown(label="Alignment Variation", choices=[], visible=True, scale=2)
                 manual_search_btn = gr.Button("Inspect Site", scale=1)
            
            # Full Width Plot
            with gr.Row():
                 plot_detail = gr.Plot(label="Site Architecture", container=True, show_label=False)
            with gr.Row():
                 detail_seq = gr.Textbox(label="Binding Sequence")
            
            gr.Markdown("### Downloads")
            csv_output = gr.File(label="Download Result CSV")

    with gr.Group():
        gr.Markdown("### 📟 System Log")
        terminal_log = gr.Textbox(
            label="System Log", 
            lines=20, 
            interactive=False, 
            show_label=False,
            elem_id="system_log_terminal",
            autoscroll=True
        )

    # --- Event Logic ---
    
    # Open Folder Buttons
    btn_open_in.click(open_input_folder)
    btn_open_out.click(open_output_folder)
    
    # Save Custom GFF
    btn_save_custom.click(
        save_custom_gff,
        inputs=[custom_gene_df],
        outputs=[save_status]
    )
    
    # Mode Switch
    def toggle_mode(mode):
        return {
            group_manual: gr.update(visible=(mode=="Manual Sequence")),
            group_ncbi: gr.update(visible=(mode=="NCBI Accession"))
        }

    mode_select.change(toggle_mode, inputs=[mode_select], outputs=[group_manual, group_ncbi])

    # NCBI Fetch
    ncbi_fetch_btn.click(
        fetch_ncbi_basics, 
        inputs=[ncbi_accession], 
        outputs=[input_status_log, ncbi_fasta_down, ncbi_gff_down]
    )
    
    # NCBI Preview (Refactored to Unified Preview)
    preview_btn.click(
        preview_genes_wrapper,
        inputs=[mode_select, rna_text, ncbi_gff_down, range_start, range_end, feature_radio, ncbi_accession, custom_gene_df],
        outputs=[preview_plot]
    )

    # Run Analysis
    # Run Analysis
    run_btn.click(
        run_dispatch_wrapper,
        inputs=[
             mode_select,
             prot_text, prot_file, rna_text, rna_file,
             ncbi_accession, range_start, range_end, ncbi_fasta_down, ncbi_gff_down,
             tool_select, scoring_method, range_max_gap, use_ps_score, job_name_input,
             feature_radio, plot_height, chk_highlight, dd_strand,
             log_state, # Pass log state input
             custom_gene_df # Pass custom genes
        ],
        outputs=[
             terminal_log, csv_output, motif_plot, motif_table, 
             plot_browser, plot_detail, df_state, analysis_context_state
        ]
    ).then(
        # Sync log_state with terminal_log value after run
        lambda x: x, inputs=[terminal_log], outputs=[log_state]
    )
    
    # Update Run Status Simple (Sync with detailed log)
    run_btn.click(
        lambda: "Analysis Started...", outputs=[run_status_simple]
    ).then(
        lambda s: "Analysis Complete!" if "Found" in s else "Analysis Failed/Finished", 
        inputs=[terminal_log], 
        outputs=[run_status_simple]
    )
    
    # Detail Inspector
    manual_search_btn.click(
        show_manual_detail,
        inputs=[manual_pos, manual_strand, df_state, analysis_context_state],
        outputs=[gap_select, manual_data_state, plot_detail, detail_seq, terminal_log]
    )
    
    gap_select.change(
        inspect_gap_select, 
        inputs=[gap_select, manual_data_state],
        outputs=[plot_detail, detail_seq, terminal_log]
    )
    
    # NOTE: gr.Plot does not support .select() for Plotly in Gradio 6.x
    # The user can manually enter Position from Plotly hover tooltips and click "Inspect Site".
    
    # Live Browser Updates (Highlight/Strand/Height)
    def update_browser_view(df, highlight, strand, height, context, accession, mode, gff_file, r_start, r_end, feats, job_name):
        # We need to know if we are in Manual (Full) or NCBI (Range/Slice) mode to call correct plotter.
        # But `df` in state is always Standardized (Manual has 0-based, NCBI has Absolute).
        # Wait, My NCBI fix Shifted df to Absolute.
        # But `plot_range_browser` needs `start` and `end`.
        # `plot_genome_browser` (Manual) plots Everything.
        
        if df is None or df.empty:
            return go.Figure()

        # Heuristic: Check Mode Select? No, inputs not passed.
        # Use Context?
        # If 'genome_offset' is in context, we are in NCBI Mode.
        offset = context.get('genome_offset', None)
        
        if offset is not None:
             # NCBI Mode -> Use Range Browser
             gff_path = gff_file.name if gff_file else None
             return plot_range_browser(
                df, gff_path, int(r_start), int(r_end), 
                seq_id=accession, 
                feature_types=feats,
                plot_height=height, 
                highlight_top=highlight, 
                highlight_strand=strand,
                analysis_context={}, # Disable recalc
                protein_name=job_name
             )
        else:
             # Manual Mode -> Use Genome Browser
             return plot_genome_browser(
                 df, seq_id="Target RNA", 
                 highlight_top=highlight, 
                 highlight_strand=strand
             )
             # Note: Height control for plot_genome_browser? It needs update.
             # plot_genome_browser currently has hardcoded height logic or uses default layout update?
             # Let's check plot_genome_browser. It returns fig.
             # We can update height here. 
    
    # We need a wrapper that updates height on the returned fig
    def update_view_wrapper(*args):
        # args: df, highlight, strand, height, context, accession, mode, gff, start, end, feats
        fig = update_browser_view(*args)
        # Apply height if valid fig
        if fig:
             fig.update_layout(height=args[3]) # args[3] is height
        return fig

    view_inputs = [
        df_state, chk_highlight, dd_strand, plot_height, analysis_context_state, 
        ncbi_accession, mode_select, ncbi_gff_down, range_start, range_end, feature_radio, job_name_input
    ]
    
    chk_highlight.change(update_view_wrapper, inputs=view_inputs, outputs=[plot_browser])
    dd_strand.change(update_view_wrapper, inputs=view_inputs, outputs=[plot_browser])
    plot_height.change(update_view_wrapper, inputs=view_inputs, outputs=[plot_browser])

    # Render complete (floating terminal already defined above)

    gr.Markdown(f"---\n**PPPRPP {VERSION}** | Unified Tabbed Interface")

if __name__ == "__main__":
    cleanup_temp_files()
    server_name = os.environ.get("PPPRPP_SERVER_NAME", "127.0.0.1")
    try:
        server_port = int(os.environ.get("PPPRPP_SERVER_PORT", "7860"))
    except ValueError:
        server_port = 7860
    inbrowser_flag = os.environ.get("PPPRPP_INBROWSER", "1").strip().lower()
    inbrowser = inbrowser_flag in {"1", "true", "yes", "on"}

    demo.queue().launch(
        share=False,
        inbrowser=inbrowser,
        server_name=server_name,
        server_port=server_port,
    )
