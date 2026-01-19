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

# --- Constants ---
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "data")
OUTPUT_DIR = os.path.join(DATA_DIR, "output")
INPUT_DIR = os.path.join(DATA_DIR, "input")

# Ensure directories exist
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(INPUT_DIR, exist_ok=True)

TMP_DIR = OUTPUT_DIR

VERSION = "v0.1.0"

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

# Import package modules
from ppr_analysis import core, io as ppr_io, visualization, config, motif_generation, ncbi

# Check Gradio Version
if int(gr.__version__.split('.')[0]) < 4:
    raise RuntimeError(f"Gradio version {gr.__version__} is too old. Please install gradio>=4.0.0")
print(f"Gradio Version: {gr.__version__}")

custom_css = """
#debug_console {
    position: fixed;
    bottom: 10px;
    right: 10px;
    width: 400px;
    height: 300px;
    z-index: 9999;
    opacity: 0.9;
    border: 1px solid #ccc;
    box-shadow: 0 0 10px rgba(0,0,0,0.2);
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
        df = pd.read_csv(path, usecols=usecols)
        return df
    except Exception as e:
        print(f"Error loading optimized CSV: {e}")
        return pd.DataFrame()

        return pd.DataFrame()

def parse_motif_file_to_records(motifs_file, seq_id="Unknown"):
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
        motif_data = {'motif_labels': motif_labels, 'ppr_scores': ppr_scores}
            
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
             motif_df, motif_fig = visualize_motifs(motif_path, analysis_title=f"{job_name} Motifs")
        
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
    sequence = core.dna_to_rna(sequence)

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
            # Use the Web API version to avoid local binary issues
            pkl_path = motif_generation.run_ps_scan_api(tmp_fasta_path, output_dir)
            seq_id = os.path.basename(pkl_path).replace('.pkl', '')
            motifs_file = Path(output_dir) / f"{seq_id}_motifs.txt"
            
        if not pkl_path:
             return "Error: Motif discovery failed.", None, None, None, None, None, None, None, None, None

    except Exception as e:
        return f"Error during discovery: {e}", None, None, None, None, None, None, None, None, None
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
        return f"Error loading motifs: {e}", None, None, None, None, None, None, None, None, None

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
        
        # Reverse
        rev_seq = core.get_reverse_complement(rna_seq)
        scores_rev = core.analyze_ppr_data(motif_data, rev_seq, min_front, min_back, max_gap, predicted_scores, ppr_name, detailed_output=False, target_analysis=target_analysis_config)
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
                
    except Exception as e:
        import traceback
        traceback.print_exc()
        return f"Error during analysis: {e}", None, None, None, None, None, None, None, None, None
    if not all_scores["Position_Relative"]:
        return "Error: No binding sites found.", None, None, None, None, None, None, None, None, None

    # Create DataFrame
    df_output = pd.DataFrame(all_scores)
    
    # Add Analysis Type
    df_output['Analysis_Type'] = analysis_type
    
    # Filter for significant hits (optional, but good for default view)
    # For now, return all or top hits. Let's return all.
    
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
    
    # Save to CSV using the new descriptive name
    # Ensure it ends with _results.csv logic for past_results glob
    # Actually, we appended _results to final_job_name, so just add .csv
    csv_path = os.path.join(TMP_DIR, f"{final_job_name}.csv")
    df_output.to_csv(csv_path, index=False)
    
    # Save accompanying motifs file for Load Result
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
        print(f"DEBUG: Motifs file missing or invalid. Source: {motifs_file}")

    # 5. Visualization (Genome Browser)
    fig = go.Figure()
    if not skip_plots:
        fig = plot_genome_browser(df_output, seq_id, analysis_type, highlight_top=highlight_top, highlight_strand=highlight_strand)

    
    # 6. Visualize Motifs
    motif_df, motif_fig = pd.DataFrame(), go.Figure()
    if motifs_file and os.path.exists(motifs_file):
        friendly_title = f"Method: {friendly_method} | Gap: {int(max_gap_size)} | Score: {'On' if use_ps_score else 'Off'}"
        motif_df, motif_fig = visualize_motifs(motifs_file, analysis_title=friendly_title)
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

def plot_genome_browser(df, seq_id="Sequence", analysis_type="Analysis", highlight_top=False, highlight_strand="Both"):
    """Helper to plot the genome browser."""
    # Try to get analysis type from DF if not provided or generic
    if (analysis_type == "Analysis" or analysis_type is None) and 'Analysis_Type' in df.columns and not df.empty:
        analysis_type = df['Analysis_Type'].iloc[0]
        
    readable_condition = parse_analysis_type(analysis_type)
    
    fig = make_subplots(
        rows=2, cols=1, 
        shared_xaxes=True, 
        vertical_spacing=0.08, 
        subplot_titles=("Forward Strand Binding", "Reverse Strand Binding")
    )
    
    unique_ppr_names = df['PPR_Name'].unique()
    colors = px.colors.qualitative.Plotly
    
    # Calculate global max score first to ensure consistent Y-axis
    max_score = df['PPR_Binding_Score'].abs().max() if not df.empty else 10
    y_limit = max_score * 1.1

    for i, name in enumerate(unique_ppr_names):
        color = colors[i % len(colors)]
        
        # Forward
        df_fwd = df[(df['PPR_Name'] == name) & (df['Strand'] == "Forward")]
        if not df_fwd.empty:
            custom_data = np.stack((df_fwd.index, df_fwd['PPR_Binding_Score']), axis=-1)
            fig.add_trace(go.Scatter(
                x=df_fwd['Genome_Position'],
                y=df_fwd['PPR_Binding_Score'],
                mode='markers',
                name=f"{name} (+)",
                marker=dict(color=color, size=6, symbol='circle'),
                customdata=custom_data, 
                hovertemplate="Pos: %{x}<br>Score: %{y:.2f}<br>Strand: +<extra></extra>"
            ), row=1, col=1)
            
        # Reverse
        df_rev = df[(df['PPR_Name'] == name) & (df['Strand'] == "Reverse")]
        if not df_rev.empty:
            custom_data = np.stack((df_rev.index, df_rev['PPR_Binding_Score']), axis=-1)
            fig.add_trace(go.Scatter(
                x=df_rev['Genome_Position'],
                y=df_rev['PPR_Binding_Score'], 
                mode='markers',
                name=f"{name} (-)",
                marker=dict(color=color, size=6, symbol='x'),
                customdata=custom_data,
                hovertemplate="Pos: %{x}<br>Score: %{y:.2f}<br>Strand: -<extra></extra>"
            ), row=2, col=1)

    # Highlight Top Score
    if highlight_top and not df.empty:
        # Filter by strand if requested
        if highlight_strand == "Forward":
            df_search = df[df['Strand'] == 'Forward']
        elif highlight_strand == "Reverse":
            df_search = df[df['Strand'] == 'Reverse']
        else:
            df_search = df
            
        if not df_search.empty:
            top_idx = df_search['PPR_Binding_Score'].abs().idxmax()
            top_row = df.loc[top_idx] # Access from original df to be safe with indices
            
            t_pos = top_row['Genome_Position']
            t_score = top_row['PPR_Binding_Score']
            t_strand = top_row['Strand']
            
            # Determine row based on strand
            row_idx = 1 if t_strand == "Forward" else 2
            
            fig.add_trace(go.Scatter(
                x=[t_pos], y=[t_score],
                mode='markers',
                marker=dict(symbol='star', size=15, color='gold', line=dict(width=2, color='black')),
                name=f"Top Score ({t_strand})",
                hoverinfo='text',
                hovertext=f"TOP SCORE<br>Pos: {t_pos}<br>Score: {t_score:.2f}"
            ), row=row_idx, col=1)

    fig.update_layout(
        title=f"Genome Browser - {seq_id}<br><sup>{readable_condition}</sup>",
        height=500,
        hovermode="closest",
        margin=dict(t=80, b=50, l=50, r=50) # Increased top margin for subtitle
    )
    
    fig.update_xaxes(title_text="Genomic Position", row=2, col=1)
    fig.update_yaxes(title_text="Fwd Score", row=1, col=1, range=[0, y_limit])
    fig.update_yaxes(title_text="Rev Score", row=2, col=1, range=[0, y_limit], autorange="reversed")
    
    return fig

def plot_range_browser(df, gff_path, start, end, seq_id="Sequence", feature_types=None, plot_height=500, highlight_top=False, highlight_strand="Both"):
    """Plots binding scores and gene annotations for a specific range."""
    # Extract analysis type from DF
    analysis_type = "Analysis"
    if 'Analysis_Type' in df.columns and not df.empty:
        analysis_type = df['Analysis_Type'].iloc[0]
    
    readable_condition = parse_analysis_type(analysis_type)
    
    if feature_types is None:
        feature_types = ["gene", "CDS", "tRNA", "rRNA", "exon"]
        
    fig = make_subplots(
        rows=3, cols=1, 
        shared_xaxes=True, 
        vertical_spacing=0.08, 
        row_heights=[0.4, 0.2, 0.4],
        subplot_titles=("Forward Strand Binding", "Gene Annotations", "Reverse Strand Binding")
    )
    
    # Calculate global max score first for consistent Y-axis
    max_score = df['PPR_Binding_Score'].abs().max() if not df.empty else 10
    y_limit = max_score * 1.1

    # 1. Plot Binding Scores (Row 1: Forward, Row 3: Reverse)
    df_range = df[(df['Genome_Position'] >= start) & (df['Genome_Position'] <= end)]
    
    unique_ppr_names = df['PPR_Name'].unique()
    colors = px.colors.qualitative.Plotly
    
    for i, name in enumerate(unique_ppr_names):
        color = colors[i % len(colors)]
        
        # Forward
        df_fwd = df_range[(df_range['PPR_Name'] == name) & (df_range['Strand'] == "Forward")]
        if not df_fwd.empty:
            fig.add_trace(go.Scatter(
                x=df_fwd['Genome_Position'],
                y=df_fwd['PPR_Binding_Score'],
                mode='markers',
                name=f"{name} (+)",
                marker=dict(color=color, size=8, symbol='circle'),
                hovertemplate="Pos: %{x}<br>Score: %{y:.2f}<br>Strand: +<extra></extra>"
            ), row=1, col=1)
            
        # Reverse
        df_rev = df_range[(df_range['PPR_Name'] == name) & (df_range['Strand'] == "Reverse")]
        if not df_rev.empty:
            fig.add_trace(go.Scatter(
                x=df_rev['Genome_Position'],
                y=df_rev['PPR_Binding_Score'], 
                mode='markers',
                name=f"{name} (-)",
                marker=dict(color=color, size=8, symbol='x'),
                hovertemplate="Pos: %{x}<br>Score: %{y:.2f}<br>Strand: -<extra></extra>"
            ), row=3, col=1)

    # Highlight Top Score in Range
    if highlight_top and not df_range.empty:
         # Filter by strand if requested
        if highlight_strand == "Forward":
            df_search = df_range[df_range['Strand'] == 'Forward']
        elif highlight_strand == "Reverse":
            df_search = df_range[df_range['Strand'] == 'Reverse']
        else:
            df_search = df_range
            
        if not df_search.empty:
            top_idx = df_search['PPR_Binding_Score'].abs().idxmax()
            top_row = df_range.loc[top_idx] 
            
            t_pos = top_row['Genome_Position']
            t_score = top_row['PPR_Binding_Score']
            t_strand = top_row['Strand']
            
            # Row 1 for Fwd, Row 3 for Rev
            row_idx = 1 if t_strand == "Forward" else 3
            
            fig.add_trace(go.Scatter(
                x=[t_pos], y=[t_score],
                mode='markers',
                marker=dict(symbol='star', size=15, color='gold', line=dict(width=2, color='black')),
                name=f"Top Score ({t_strand})",
                hoverinfo='text',
                hovertext=f"TOP SCORE<br>Pos: {t_pos}<br>Score: {t_score:.2f}"
            ), row=row_idx, col=1)

    # 2. Plot Gene Annotations from GFF3 (Row 2)
    fig.add_trace(go.Scatter(x=[start, end], y=[0, 0], mode='markers', marker=dict(opacity=0), showlegend=False, hoverinfo='skip'), row=2, col=1)
    
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
                type="path", path=path, line=dict(color=color, width=1), fillcolor=color, opacity=0.6, row=2, col=1
            )
            
            # Add label
            fig.add_annotation(
                x=(f_start + f_end) / 2,
                y=y_pos,
                text=f_name,
                showarrow=False,
                font=dict(color="black", size=10),
                yshift=15 if f_strand == '+' else -15,
                row=2, col=1
            )
            
    # Add title with conditions
    fig.update_layout(
        title=f"{seq_id}<br><sup>{readable_condition}</sup>",
        height=plot_height,
        hovermode="closest",
        margin=dict(t=80, b=50, l=50, r=50)
    )
            
    fig.update_yaxes(title_text="Rev Score", row=3, col=1, range=[0, y_limit], autorange="reversed")
    
    return fig


def visualize_motifs(motifs_file, analysis_title=None):
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
        
        # Layout settings
        title_text = f"Discovered Motifs (PPRFinder) - {fasta_name}"
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
        
        for i, label in enumerate(motif_labels):
            # X-axis: Nucleotide + Index + Motif Number
            nuc = segment[i] if i < len(segment) else "-"
            
            # Construct label
            # If gap, explicitly show GAP
            if label == '--':
                m_num_str = "(Gap)"
                text_label = "GAP"
                color = color_map['--']
            else:
                m_num = motif_numbers[i]
                m_num_str = f"#{m_num}" if m_num else ""
                text_label = f"{m_num_str}<br>{label}"
                color = color_map.get(label, '#bcbd22')

            # X label: Position (1-based relative) <br> Nucleotide <br> Motif #
            x_labels.append(f"{i+1}<br><b>{nuc}</b><br><span style='font-size:10px'>{m_num_str}</span>")
            
            # Colors list
            colors.append(color)
            
            # Text on bar (Motif Type)
            text_labels.append(label if label != '--' else 'GAP')

        fig = go.Figure()
        
        # Bar chart for Binding Scores
        # If score is 0, bar won't show. Add a small value for GAP visibility?
        # Or just rely on x-axis labels.
        # Let's add a separate trace for Gaps if needed, but text_labels help.
        
        fig.add_trace(go.Bar(
            x=x_labels, 
            y=scores, 
            text=text_labels,
            textposition='auto', 
            marker_color=colors,
            name="Binding Score"
        ))
        
        # Line chart for Confidence Scores (if available)
        if confidence_scores and len(confidence_scores) == len(scores):
            fig.add_trace(go.Scatter(
                x=x_labels,
                y=confidence_scores,
                mode='lines+markers',
                name='Motif Confidence',
                line=dict(color='black', width=2, dash='dot'),
                yaxis='y2'
            ))
        
        # Get Analysis Type if available
        analysis_type_str = ""
        # Improved formatting for Analysis Type and Site details
        if 'Analysis_Type' in row and pd.notna(row['Analysis_Type']):
             atype = row['Analysis_Type']
             # Parse: method_gap_ps
             parts = atype.split('_')
             
             # Method
             method_map = {'kobayashi': 'PLS-type', 'yan': 'P-type'}
             method_name = method_map.get(parts[0], parts[0].capitalize())
             
             # Gap
             gap_str_raw = parts[1] if len(parts) > 1 else ""
             # Check ps
             # User requested clearer label for Motif Confidence Score (+PS)
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
        x_val = evt.value['x']
        y_val = evt.value['y']
        
        # Approximate matching for floats
        mask = (df['Genome_Position'] == x_val) & \
               ((np.isclose(df['PPR_Binding_Score'], y_val)) | (np.isclose(df['PPR_Binding_Score'], -y_val)))
        
        row = df[mask]
        
        if row.empty:
             return gr.update(choices=[]), gr.update(value=None), go.Figure(layout=dict(title="No data found for clicked point")), "", f"No data found for x={x_val}, y={y_val}"
            
        best_row = row.iloc[0]
        # Strand is capitalized in DF: "Forward", "Reverse"
        strand = best_row['Strand']
        
        # Delegate to show_manual_detail to handle Gap Selection logic (Dynamic Mode if context available)
        return show_manual_detail(str(x_val), strand, df, analysis_context)
        
    except Exception as e:
        print(f"Error in click handler: {e}")
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
            
        target_pos = int(manual_pos)
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
        
        # Handle strand
        if manual_strand == "Reverse":
             seq_to_analyze = core.get_reverse_complement(rna_seq)
             rev_pos = len(rna_seq) - 1 - target_pos
             calc_pos = rev_pos
        else:
             seq_to_analyze = rna_seq
             calc_pos = target_pos
             
        logs.append(f"Analyzing single site. Calc Pos: {calc_pos}, Max Gap: {max_gap}")
             
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
                 return gr.update(choices=[], value="Error: Detailed data missing and Context unavailable."), None, go.Figure(), "", "\n".join(logs)
            else:
                 # This branch shouldn't be reached if logic above (if analysis_context: use_dynamic=True) is correct.
                 # Wait, logic above checks: if 'motif_data' in analysis_context...
                 # If analysis_context exists but missing keys?
                 logs.append("Context exists but might be malformed or `use_dynamic` logic failed.")
                 # Force return error or try again?
                 return gr.update(choices=[], value="Error: Lite CSV loaded but Dynamic Analysis failed."), None, go.Figure(), "", "\n".join(logs)
                 
        target_data = subset.to_dict('records')
        logs.append(f"Found {len(target_data)} matching rows in DF.")

    if not target_data:
        logs.append("No target data found.")
        return gr.update(choices=[], value="No data"), None, go.Figure(), "", "\n".join(logs)

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
    # Outputs: [status_output, ncbi_fasta_download, ncbi_bedgraph_download, ncbi_bed_download, ncbi_gff_download, ncbi_download, ncbi_log_output]
    def make_output(status, fasta=None, bedgraph=None, bed=None, gff=None, csv=None):
        return (
            status, 
            fasta, 
            bedgraph,
            bed,
            gff,
            csv,
            "\n".join(log_buffer)
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
        result_msg, csv_path, motif_fig, _, fig_browser, _, df, _, _, _ = run_unified_analysis(
            protein_text, protein_file, rna_text, None, tool_choice, scoring_method, max_gap_size, use_ps_score, job_name, skip_plots=True
        )
        log(f"Analysis result: {result_msg}")
        
        # Optimize for UI: Truncate DataFrame if too large
        ui_df = df.head(1000) if df is not None else None
        if df is not None and len(df) > 1000:
            log(f"Truncating DataFrame for UI display (showing 1000 of {len(df)} rows). Full data in CSV.")
            
        yield make_output("Analysis Finished.", fasta=fasta_path, csv=csv_path)
        
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

# --- UI Definition ---

with gr.Blocks(title="P-type PPR Prediction Program (PPPRPP)") as demo:
    gr.Markdown(f"# P-type PPR Prediction Program (PPPRPP) <span style='font-size: 0.5em; font-weight: normal; color: #777;'>{VERSION}</span>")
    gr.Markdown("Input your Protein sequence (to discover motifs) and RNA sequence (target) to analyze binding sites.")
    
    # Debug Console (Defined early for variable access)
    # Debug Console (Defined early for variable access, but rendered later)
    with gr.Group(elem_id="debug_console", render=False) as debug_group:
        gr.Markdown("### Debug Log")
        debug_console = gr.Textbox(label="Logs", lines=12, max_lines=20, interactive=False, elem_id="debug_log_text")
    
    with gr.Row():
        with gr.Column():
            gr.Markdown("### 1. Protein Input (Motif Discovery)")
            prot_text = gr.Textbox(label="Protein FASTA", lines=5, placeholder=f"Default: >PPR791(Rf1a)_IR24\\n{DEFAULT_PROTEIN[:50]}...")
            prot_file = gr.File(label="Or Upload Protein File")
            tool_select = gr.Dropdown(choices=["PPRFinder", "PS_Scan"], value="PPRFinder", label="Discovery Tool")
            
        with gr.Column():
            gr.Markdown("### 2. RNA Input (Target Analysis)")
            rna_text = gr.Textbox(label="RNA FASTA", lines=5, placeholder=f"Default: >atp6-orf79_BT\\n{DEFAULT_RNA[:50]}...")
            rna_file = gr.File(label="Or Upload RNA File")
            
    with gr.Row():
        with gr.Column():
            gr.Markdown("### 3. Settings & Run")
            scoring_method = gr.Radio(["P-type", "PLS-type"], value="P-type", label="Scoring Method")
            max_gap_size = gr.Number(value=0, label="Gap Size (0 = No Gaps)")
            use_ps_score = gr.Checkbox(label="Use Motif Confidence Score (+Score)", value=True)
            job_name_input = gr.Textbox(value="my_analysis", label="Job Name")
            
            # Result Management
            past_results = list_past_results()
            with gr.Row():
                result_selector = gr.Dropdown(label="Past Results", choices=past_results, value=None, interactive=True, scale=4)
                load_btn = gr.Button("Load Result", variant="secondary", scale=1)
            
            run_btn = gr.Button("Run Unified Analysis", variant="primary")
            
        with gr.Column():
            status_output = gr.Textbox(label="Status", lines=15)
            
    # Outputs
    with gr.Row():
        with gr.Column():
            gr.Markdown("### Motif Discovery Results")
            motif_plot = gr.Plot(label="Motif Architecture")
            motif_table = gr.DataFrame(label="Motif Details")
            
    with gr.Row():
        chk_highlight_main = gr.Checkbox(label="Highlight Top Score", value=False)
        dd_strand_main = gr.Dropdown(choices=["Forward", "Reverse", "Both"], value="Both", label="Highlight Strand", interactive=True)
        
    with gr.Row():
        plot_browser = gr.Plot(label="Genome Browser (Click point to see details)")
    
    with gr.Row():
        with gr.Column(scale=1):
            gr.Markdown("### Detail Inspector")
            with gr.Row():
                manual_pos = gr.Number(label="Genome Position", precision=0)
                manual_strand = gr.Radio(["Forward", "Reverse"], value="Forward", label="Strand")
            
            # Dropdown for Gap Selection (initially hidden/empty)
            gap_select = gr.Dropdown(label="Select Alignment Variation", choices=[], visible=True)
            
            manual_btn = gr.Button("Search Site")
        with gr.Column(scale=3):
            plot_detail = gr.Plot(label="Selected Site Detail")
            detail_seq = gr.Textbox(label="Binding Sequence", interactive=True)
        
    gr.Markdown("### Download Results")
    with gr.Row():
        csv_output = gr.File(label="Download CSV")
        
    # State
    # State
    df_state = gr.State()
    analysis_context_state = gr.State({})
    manual_data_state = gr.State({}) # Store detail data for gap selection
    ncbi_data_state = gr.State({}) # For NCBI tab Gap Selection
    
    with gr.Tab("NCBI Analysis"):
        gr.Markdown("### Import Target Sequence from NCBI and Analyze")
        gr.Markdown("Enter an NCBI Accession ID to fetch the sequence and analyze binding sites. Download the generated FASTA and BedGraph files to view them in the [IGV Desktop Application](https://igv.org/doc/desktop/).")
        
        with gr.Row():
            ncbi_accession = gr.Textbox(label="NCBI Accession ID", placeholder=f"Default: {DEFAULT_ACCESSION}", scale=2)
            chk_generate_bed = gr.Checkbox(label="Generate IGV Files", value=False, scale=1) 
            ncbi_btn = gr.Button("Fetch & Analyze", variant="primary", scale=1)
        
        with gr.Row():
            ncbi_fasta_download = gr.File(label="Download FASTA (Reference Genome)")
            ncbi_bedgraph_download = gr.File(label="Download BedGraph (Binding Scores)")
            ncbi_bed_download = gr.File(label="Download BED (Binding Sites)") # Added
            ncbi_gff_download = gr.File(label="Download GFF3 (Gene Annotations)")
        
        with gr.Row():
            ncbi_log_output = gr.Textbox(label="Analysis Log", lines=10, interactive=False)
        
        with gr.Row():
            ncbi_download = gr.File(label="Download CSV")
            
        with gr.Row():
            with gr.Row():
                ncbi_result_selector = gr.Dropdown(label="Load Past Result (NCBI)", choices=[], interactive=True, scale=4)
                ncbi_load_btn = gr.Button("Load Result", variant="secondary", scale=1)
            
            # Note: I need to update the event binding at the end of file to use `ncbi_result_selector` instead of `result_selector`
            # or add a NEW event for `ncbi_result_selector`.
            # The user requested "NCBI Viewer can also select past results".
            # So `ncbi_result_selector` should populate `ncbi_download` (CSV).

            
        gr.Markdown("### Range Inspector")
        gr.Markdown("Inspect a specific genomic range with gene annotations.")
        with gr.Row():
            range_start = gr.Number(label="Start Position", value=133162, precision=0)
            range_end = gr.Number(label="End Position", value=134622, precision=0)
            range_btn = gr.Button("Inspect Range")
        
        with gr.Row():
            highlight_ncbi = gr.Checkbox(label="Highlight Top Score", value=False)
            strand_ncbi = gr.Dropdown(choices=["Forward", "Reverse", "Both"], value="Both", label="Highlight Strand", interactive=True)
            
        with gr.Row():
            feature_types = gr.CheckboxGroup(
                ["gene", "CDS", "tRNA", "rRNA", "exon"], 
                value=["gene", "tRNA", "rRNA"], 
                label="Show Features"
            )
            plot_height = gr.Slider(minimum=400, maximum=1200, value=500, step=50, label="Plot Height (px)")
            
        range_plot = gr.Plot(label="Range Browser")
        # range_files = gr.File(label="Download Range Results (CSV & Interactive Plot)") # Removed
            
        def on_range_click(accession, start, end, csv_file, gff_file, features, height, highlight_top, highlight_strand):
            if not csv_file:
                return go.Figure(layout=dict(title="No CSV file found. Run analysis first."))
            try:
                # Optimized load
                cols = ['Genome_Position', 'PPR_Binding_Score', 'PPR_Name', 'Strand', 'Analysis_Type']
                df = load_csv_optimized(csv_file.name, usecols=cols)
                
                # GFF file object or path
                gff_path = gff_file.name if gff_file else None
                
                # Generate Figure
                fig = plot_range_browser(df, gff_path, int(start), int(end), accession, features, int(height), highlight_top=highlight_top, highlight_strand=highlight_strand)
                
                return fig
                
            except Exception as e:
                import traceback
                traceback.print_exc()
                return go.Figure(layout=dict(title=f"Error: {e}"))

        range_btn.click(
            on_range_click,
            inputs=[ncbi_accession, range_start, range_end, ncbi_download, ncbi_gff_download, feature_types, plot_height, highlight_ncbi, strand_ncbi],
            outputs=[range_plot]
        )
            
        gr.Markdown("### Inspect Specific Site")
        gr.Markdown("View detailed binding scores for a specific position.")
        with gr.Row():
            with gr.Column(scale=1):
                ncbi_manual_pos = gr.Number(label="Genome Position", precision=0)
                ncbi_manual_strand = gr.Radio(["Forward", "Reverse"], value="Forward", label="Strand")
                ncbi_gap_select = gr.Dropdown(label="Alignment Variation", choices=[], visible=True)
                ncbi_manual_btn = gr.Button("Show Detail")
            with gr.Column(scale=3):
                ncbi_detail_plot = gr.Plot(label="Site Detail")
                ncbi_detail_seq = gr.Textbox(label="Binding Sequence", interactive=True)
            
        def on_ncbi_detail_click(csv_file, fasta_file, pos, strand):
            if not csv_file:
                return go.Figure(layout=dict(title="No CSV file found. Run analysis first.")), "", "Error: No CSV file."
            try:
                # Optimized Loading: Find ONLY the target row
                df = load_single_site_optimized(csv_file.name, pos, strand)
                
                if df.empty:
                    return go.Figure(layout=dict(title="Site not found in results.")), "", "Site not found."
                
                # Reconstruct context for Dynamic Mode (Alignment Variation)
                fasta_path = fasta_file.name if fasta_file else None
                context = reconstruct_context_from_csv(csv_file.name, fasta_path=fasta_path)
                
                # show_manual_detail returns: dropdown_update, data_map, fig, seq, log
                gap_upd, data_st, fig, seq, log_txt = show_manual_detail(pos, strand, df, context)
                return gap_upd, data_st, fig, seq, log_txt
            except Exception as e:
                import traceback
                traceback.print_exc()
                # Return empty updates for all
                return gr.update(), {}, go.Figure(layout=dict(title=f"Error: {e}")), "", f"Error in ncbi detail: {e}"

        ncbi_manual_btn.click(
            on_ncbi_detail_click,
            inputs=[ncbi_download, ncbi_fasta_download, ncbi_manual_pos, ncbi_manual_strand],
            outputs=[ncbi_gap_select, ncbi_data_state, ncbi_detail_plot, ncbi_detail_seq, debug_console]
        )
        
    # Events
    # Events 

    run_btn.click(
        run_unified_analysis,
        inputs=[
            prot_text, prot_file,
            rna_text, rna_file,
            tool_select, scoring_method, max_gap_size, use_ps_score, job_name_input,
            chk_highlight_main, dd_strand_main
        ],
        outputs=[
            status_output, csv_output, motif_plot, motif_table, plot_browser, range_plot, df_state, analysis_context_state,
            result_selector, ncbi_result_selector
        ]
    )

    # 4. Detail Inspector logic
    manual_btn.click(
        show_manual_detail,
        inputs=[manual_pos, manual_strand, df_state, analysis_context_state],
        outputs=[gap_select, manual_data_state, plot_detail, detail_seq, debug_console] 
    )
    
    # Gap Dropdown Change Event
    # Gap Dropdown Change Event (Main Tab)
    gap_select.change(
        inspect_gap_select, 
        inputs=[gap_select, manual_data_state],
        outputs=[plot_detail, detail_seq, debug_console] # inspect_gap_select returns fig, seq, log
    )
    
    # Gap Dropdown Change Event (NCBI Tab)
    ncbi_gap_select.change(
        inspect_gap_select, 
        inputs=[ncbi_gap_select, ncbi_data_state],
        outputs=[ncbi_detail_plot, ncbi_detail_seq, debug_console]
    )
    
    # 5. Load Results logic (Main Tab)
    

    # NCBI Range Browser Logic
    range_btn.click(
        on_range_click, 
        inputs=[ncbi_accession, range_start, range_end, ncbi_download, ncbi_gff_download, feature_types, plot_height, highlight_ncbi, strand_ncbi], 
        outputs=[range_plot]
    )

    # Highlight Live Updates
    chk_highlight_main.change(
        update_main_browser_plot,
        inputs=[df_state, chk_highlight_main, dd_strand_main],
        outputs=[plot_browser]
    )
    dd_strand_main.change(
        update_main_browser_plot,
        inputs=[df_state, chk_highlight_main, dd_strand_main],
        outputs=[plot_browser]
    )
    
    highlight_ncbi.change(
        on_range_click,
        inputs=[ncbi_accession, range_start, range_end, ncbi_download, ncbi_gff_download, feature_types, plot_height, highlight_ncbi, strand_ncbi],
        outputs=[range_plot]
    )
    strand_ncbi.change(
        on_range_click,
        inputs=[ncbi_accession, range_start, range_end, ncbi_download, ncbi_gff_download, feature_types, plot_height, highlight_ncbi, strand_ncbi],
        outputs=[range_plot]
    )

    # 5. Load Results logic (Main Tab & NCBI Tab)
    
    # Common Output List
    load_outputs = [status_output, csv_output, plot_browser, motif_plot, motif_table, range_plot, df_state, analysis_context_state]
    
    # Main Tab
    result_selector.select(load_result, inputs=[result_selector], outputs=load_outputs)
    load_btn.click(load_result, inputs=[result_selector], outputs=load_outputs)
    
    # NCBI Tab
    # Use load_result to restore full analysis state even from NCBI tab
    ncbi_result_selector.select(load_result, inputs=[ncbi_result_selector], outputs=load_outputs)
    ncbi_load_btn.click(load_result, inputs=[ncbi_result_selector], outputs=load_outputs)
    
    # Populate dropdowns on load
    demo.load(list_past_results, outputs=[result_selector, ncbi_result_selector])
    
    # Reset context on load to avoid stale state
    # demo.load -> can't reset state easily without inputs.
    
    # NCBI Event
    
    # Reposition Start button row? No, just add input.
    # Where is ncbi_btn defined? Higher up.
    # Let's insert the checkbox near the button or implicitly adding it to inputs list.
    
    ncbi_btn.click(
        run_ncbi_analysis,
        inputs=[ncbi_accession, prot_text, prot_file, tool_select, scoring_method, max_gap_size, use_ps_score, job_name_input, chk_generate_bed],
        outputs=[ncbi_log_output, ncbi_fasta_download, ncbi_bedgraph_download, ncbi_bed_download, ncbi_gff_download, ncbi_download, ncbi_log_output]
    )
    
    try:
        plot_browser.select(
            update_detail_plot, 
            inputs=[df_state, analysis_context_state], 
            outputs=[gap_select, manual_data_state, plot_detail, detail_seq, debug_console]
        )
    except AttributeError:
        print("WARNING: plot_browser.select() failed. Interactive features disabled.")
    except Exception as e:
        print(f"WARNING: Failed to attach select event: {e}")
        
    # Debug Console moved to top
    # Debug Console Rendered at Bottom
    debug_group.render()


    gr.Markdown(f"---\n**PPPRPP {VERSION}** | Developed by Keisuke Igarashi | [GitHub Repository](https://github.com/igarashi309/PPPRPP_webUI)")

    
    gr.Markdown("""
    ### About this Software
    
    This software identifies PPR motifs using PPRFinder and PS_Scan (ScanProsite). Target prediction was implemented based on the frameworks described by Kobayashi et al. (PLS-type) and Yan et al. (P-type). Please cite the references below when using this software.

    ### References
    
    **PLS-type**: Kobayashi, T., Yagi, Y., & Nakamura, T. (2019). Comprehensive Prediction of Target RNA Editing Sites for PLS-Class PPR Proteins in Arabidopsis thaliana. *Plant and Cell Physiology*, 60(4), 862–874. [https://doi.org/10.1093/pcp/pcy251](https://doi.org/10.1093/pcp/pcy251)
    
    **P-type**: Yan, J., Yao, Y., Hong, S., Yang, Y., Shen, C., Zhang, Q., Zhang, D., Zou, T., & Yin, P. (2019). Delineation of pentatricopeptide repeat codes for target RNA prediction. *Nucleic Acids Research*, 47(7), 3728–3738. [https://doi.org/10.1093/nar/gkz075](https://doi.org/10.1093/nar/gkz075)
    
    **PPRFinder**: Gutmann, B., Royan, S., Schallenberg-Rüdinger, M., Lenz, H., Castleden, I. R., McDowell, R., Vacher, M. A., Tonti-Filippini, J., Bond, C. S., Knoop, V., & Small, I. D. (2020). The Expansion and Diversification of Pentatricopeptide Repeat RNA-Editing Factors in Plants. *Molecular Plant*, 13(2), 215–230. [https://doi.org/10.1016/j.molp.2019.11.002](https://doi.org/10.1016/j.molp.2019.11.002)
    
    **PS_Scan**: de Castro, E., Sigrist, C. J. A., Gattiker, A., Bulliard, V., Langendijk-Genevaux, P. S., Gasteiger, E., Bairoch, A., & Hulo, N. (2006). ScanProsite: detection of PROSITE signature matches and ProRule-associated functional and structural residues in proteins. *Nucleic Acids Research*, 34(suppl_2), W362–W365. [https://doi.org/10.1093/nar/gkl124](https://doi.org/10.1093/nar/gkl124)
    """)

if __name__ == "__main__":
    cleanup_temp_files()
    # Enable queue for long-running tasks and launch
    # Enable queue for long-running tasks and launch
    # Pass custom_css here
    demo.queue().launch(share=False, inbrowser=True, css=custom_css)
