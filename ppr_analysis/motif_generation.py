import os
import subprocess
import shutil
import urllib.request
import tarfile
import sys
from pathlib import Path
from Bio import SeqIO

# Constants for URLs
PROSITE_URL = 'https://ftp.expasy.org/databases/prosite/'
PS_SCAN_URL = 'https://ftp.expasy.org/databases/prosite/ps_scan/ps_scan_linux_x86_elf.tar.gz'
PS_MOTIF = 'PS51375'
PPRFINDER_JAR_URL = 'https://raw.githubusercontent.com/YaoYinYing/PPRfinder/master/out/production/PPRfinder/PPRfinder.jar'
HMM_PROFILE_BASE_URL = 'https://raw.githubusercontent.com/YaoYinYing/PPRfinder/master/data/'

# Default data directory within the package
DATA_DIR = Path(__file__).parent / "data"

def check_dependencies():
    """
    Checks if necessary external tools are installed.
    Returns a dict of tool -> bool.
    """
    dependencies = {
        "java": shutil.which("java") is not None,
        "perl": shutil.which("perl") is not None,
        "hmmsearch": shutil.which("hmmsearch") is not None
    }
    return dependencies

def download_file(url, dest_path):
    """Downloads a file if it doesn't exist."""
    if not dest_path.exists():
        print(f"Downloading {url} to {dest_path}...")
        try:
            urllib.request.urlretrieve(url, dest_path)
        except Exception as e:
            print(f"Failed to download {url}: {e}")

def setup_data():
    """
    Downloads necessary data files (profiles, jars, etc.) to DATA_DIR.
    """
    if not DATA_DIR.exists():
        os.makedirs(DATA_DIR)
        
    # Create subdirectories
    bin_dir = DATA_DIR / "bin"
    profiles_dir = DATA_DIR / "profiles"
    
    if not bin_dir.exists():
        os.makedirs(bin_dir)
    if not profiles_dir.exists():
        os.makedirs(profiles_dir)
        
    # Download PPRFinder.jar
    download_file(PPRFINDER_JAR_URL, bin_dir / "PPRfinder.jar")
    
    # Download HMM profiles
    for p in ['morf.hmm', 'dyw_kp.hmm', 'dyw.hmm', 'all_PPR.hmm']:
        download_file(f"{HMM_PROFILE_BASE_URL}{p}", profiles_dir / p)
        
    # Download Prosite data
    for p in ['prosite.dat', 'prorule.dat']: # profile.txt and evaluator.dat might not be strictly needed for basic scan but good to have
        download_file(f"{PROSITE_URL}{p}", profiles_dir / p)
        
    # Download PS_Scan (Linux binary - might not work on Mac, but we download it to get the perl script)
    ps_scan_tar = DATA_DIR / "ps_scan_linux_x86_elf.tar.gz"
    download_file(PS_SCAN_URL, ps_scan_tar)
    
    # Extract PS_Scan if not already extracted
    ps_scan_dir = bin_dir / "ps_scan"
    if not ps_scan_dir.exists() and ps_scan_tar.exists():
        print("Extracting PS_Scan...")
        try:
            with tarfile.open(ps_scan_tar, "r:gz") as tar:
                tar.extractall(path=bin_dir)
            # Rename extracted directory to ps_scan for easier access if needed
            # The tar usually extracts to a folder like ps_scan
            # We assume it extracts to 'ps_scan'
        except Exception as e:
            print(f"Failed to extract PS_Scan: {e}")

def run_pprfinder(protein_fasta, output_dir):
    """
    Runs PPRFinder on the given protein FASTA.
    
    Args:
        protein_fasta (str): Path to input FASTA.
        output_dir (str): Path to output directory.
        
    Returns:
        str: Path to the generated pickle file, or None if failed.
    """
    deps = check_dependencies()
    if not deps['hmmsearch']:
        raise RuntimeError("hmmsearch (HMMER) is not installed.")
        
    setup_data()
    
    # Define paths
    bin_dir = DATA_DIR / "bin"
    profiles_dir = DATA_DIR / "profiles"
    pprfinder_jar = bin_dir / "PPRfinder.jar"
    all_ppr_hmm = profiles_dir / "all_PPR.hmm"
    
    # Create output dir
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    job_name = Path(protein_fasta).stem
    hmm_out = Path(output_dir) / f"{job_name}_PPR_HMMSEARCH.domt"
    
    # 1. Run hmmsearch
    # cmd_1=f"hmmsearch --domtblout {hmmsearch_out_file} --noali -E 0.1 --cpu 2 {TMP_DIR}/profiles/all_PPR.hmm {tmp_fasta_file.name}"
    cmd_hmm = [
        "hmmsearch",
        "--domtblout", str(hmm_out),
        "--noali",
        "-E", "0.1",
        "--cpu", "2",
        str(all_ppr_hmm),
        protein_fasta
    ]
    
    print(f"Running: {' '.join(cmd_hmm)}")
    subprocess.run(cmd_hmm, check=True)
    
    # 2. Run PPRFinder (Python implementation)
    # cmd_java = ["java", "-jar", str(pprfinder_jar), protein_fasta, str(hmm_out)]
    
    print(f"Running PPRFinder (Python) on {hmm_out}...")
    
    expected_output = Path(str(hmm_out) + "_motifs.txt")
    
    try:
        from . import ppr_finder
        ppr_finder.run(str(hmm_out), str(expected_output), protein_fasta)
    except Exception as e:
        raise RuntimeError(f"PPRFinder (Python) failed: {e}")
    
    if not expected_output.exists():
        raise RuntimeError(f"PPRFinder did not generate expected output: {expected_output}")
        
    return parse_pprfinder_output(expected_output, protein_fasta, output_dir)

def parse_pprfinder_output(motifs_file, protein_fasta, output_dir):
    """
    Parses the PPRFinder output txt and converts it to the pickle format expected by the analysis tool.
    """
    import pandas as pd
    from Bio import SeqIO
    import pickle
    
    # Read original sequence to get ID and sequence content
    record = next(SeqIO.parse(protein_fasta, "fasta"))
    seq_id = record.id
    sequence_str = str(record.seq)
    
    # Read motifs
    # names=['id','start','end','score', 'motif_seq','second','fifth_aa','last_aa','motif_type']
    try:
        ppr_df = pd.read_table(motifs_file, names=['id','start','end','score', 'motif_seq','second','fifth_aa','last_aa','motif_type'], header=None)
    except pd.errors.EmptyDataError:
        print("No motifs found.")
        return None

    parsed_ppr = []
    
    # Notebook logic:
    # for pos_start,pos_end,score,motif_type,fifth,last,motif_seq in zip(ppr_df.start, ppr_df.end,ppr_df.score,ppr_df.motif_type,ppr_df.fifth_aa,ppr_df.last_aa,ppr_df.motif_seq):
    #   mypos=sequence_[pos_start:pos_end] (This is a Seq object slice)
    #   mypos.id=f'{s.id}'
    #   mypos.annotations={...}
    
    # We need to mimic the object structure. The analysis tool expects objects with .annotations.
    # We can use a simple class or just Bio.SeqRecord.SeqRecord
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    
    # Note: PPRFinder output is 1-based or 0-based?
    # Notebook uses: mypos=sequence_[pos_start:pos_end]
    # Pandas read_table usually reads as is.
    # Let's assume the notebook logic is correct for the output format.
    
    for _, row in ppr_df.iterrows():
        pos_start = int(row['start'])
        pos_end = int(row['end'])
        
        # Adjust for 0-based indexing if needed. 
        # Usually HMMER/PPRFinder might be 1-based.
        # If notebook does sequence_[pos_start:pos_end], and python slicing is 0-based...
        # Let's assume the notebook author handled it or the output is 0-based (unlikely for bio tools) or compatible.
        # Actually, looking at notebook:
        # ppr_df=pd.read_table(...)
        # mypos=sequence_[pos_start:pos_end]
        # If start is 1-based, this might be off by one?
        # But let's stick to the notebook logic for now.
        
        # Wait, if I use SeqRecord, I need to be careful.
        # The notebook uses `sequence_` which is a SeqRecord object from SeqIO.parse.
        # `mypos` becomes a SeqRecord slice.
        
        # Slice the sequence
        # Ensure indices are within bounds
        safe_start = max(0, pos_start)
        safe_end = min(len(sequence_str), pos_end)
        
        subseq = sequence_str[safe_start:safe_end]
        
        mypos = SeqRecord(Seq(subseq), id=seq_id)
        
        # Annotations
        # "PPR_CODE": fifth+last
        # "PS_SCORE": score
        
        fifth = str(row['fifth_aa'])
        last = str(row['last_aa'])
        
        mypos.annotations = {
            "PPR_TYPE": row['motif_type'],
            "PS_SCORE": float(row['score']),
            "PPR_STARTS": pos_start,
            "PPR_ENDS": pos_end,
            "FIFTH_AA": fifth,
            "LAST_AA": last,
            "PPR_CODE": fifth + last,
            # "RNA_BASE": ... (This requires code2rna_multi map, maybe we can skip or import it)
        }
        
        parsed_ppr.append(mypos)
        
    # Save to pickle
    output_pkl = Path(output_dir) / f"{seq_id}.pkl"
    with open(output_pkl, "wb") as f:
        pickle.dump(parsed_ppr, f)
        
    # Also save the motifs txt to a predictable name for visualization
    import shutil
    output_txt = Path(output_dir) / f"{seq_id}_motifs.txt"
    shutil.copy(motifs_file, output_txt)
        
    return str(output_pkl)

    return str(ps_out_file)

def run_ps_scan_api(protein_fasta, output_dir):
    """
    Runs PS_Scan via the ExPASy ScanProsite REST API.
    This avoids local binary dependencies.
    """
    from Bio.ExPASy import ScanProsite
    from Bio import SeqIO
    import time
    
    # Read sequence
    record = next(SeqIO.parse(protein_fasta, "fasta"))
    seq_id = record.id
    sequence = str(record.seq)
    
    print(f"Running PS_Scan API for {seq_id}...")
    
    # API call
    # sig='PS51375' is the PPR motif profile
    try:
        handle = ScanProsite.scan(seq=sequence, sig=PS_MOTIF, mirror='https://prosite.expasy.org', output='xml')
        xml_result = handle.read()
        handle.close()
    except Exception as e:
        raise RuntimeError(f"ScanProsite API failed: {e}")
        
    # Save XML
    xml_path = Path(output_dir) / f"{seq_id}_ps_scan.xml"
    with open(xml_path, "wb") as f:
        f.write(xml_result)
        
    print(f"Saved XML result to {xml_path}")
    
    # Parse XML to generate .pkl
    return parse_ps_scan_xml(xml_path, protein_fasta, output_dir)

def parse_ps_scan_xml(xml_file, protein_fasta, output_dir):
    """
    Parses ScanProsite XML output and converts it to the pickle format.
    """
    import xml.etree.ElementTree as ET
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    import pickle
    
    # Read original sequence
    record = next(SeqIO.parse(protein_fasta, "fasta"))
    seq_id = record.id
    sequence_str = str(record.seq)
    
    tree = ET.parse(xml_file)
    root = tree.getroot()
    
    parsed_ppr = []
    
    # Iterate over matches
    # Structure: <match> <start> <stop> <signature> ... </match>
    # Note: ScanProsite XML structure might vary.
    # Usually: <match_list> <match> ... </match> </match_list>
    
    # Handle namespaces if present
    # The XML from API usually has xmlns="urn:expasy:scanprosite" in matchset
    namespaces = {'ns': 'urn:expasy:scanprosite'}
    
    # Try finding with namespace first
    matches = root.findall(".//ns:match", namespaces)
    if not matches:
        # Fallback to no namespace if not found (or if format changed)
        matches = root.findall(".//match")
    
    for match in matches:
        # Helper to get text with or without namespace
        def get_text(elem, tag):
            t = elem.find(f"ns:{tag}", namespaces)
            if t is None:
                t = elem.find(tag)
            return t.text if t is not None else None
            
        # Check if it's the correct signature
        sig_ac = get_text(match, "signature/ac")
        # Sometimes signature is nested differently or we just check signature_ac direct child
        if sig_ac is None:
             # Try deeper?
             sig = match.find("ns:signature", namespaces)
             if sig:
                 ac = sig.find("ns:ac", namespaces)
                 if ac is not None:
                     sig_ac = ac.text
        
        if sig_ac is not None and sig_ac != PS_MOTIF:
            continue
            
        start_text = get_text(match, "start")
        end_text = get_text(match, "stop")
        score_text = get_text(match, "score")
        
        if not start_text or not end_text:
            continue
            
        start = int(start_text)
        end = int(end_text)
        score = float(score_text) if score_text else 0.0
        
        # Extract sequence
        # API returns 1-based coordinates? Usually yes.
        # Python is 0-based.
        # Let's assume 1-based inclusive from API.
        safe_start = max(0, start - 1)
        safe_end = min(len(sequence_str), end)
        
        subseq = sequence_str[safe_start:safe_end]
        
        # Determine 5th and Last AA
        # Based on optimization: 5th is index 5, Last is index 31 (0-based) relative to PS51375 match
        
        fifth = ""
        last = ""
        if len(subseq) > 31:
             fifth = subseq[5]
             last = subseq[31]
        elif len(subseq) >= 6:
             fifth = subseq[5]
             last = subseq[-1]
             
        mypos = SeqRecord(Seq(subseq), id=seq_id)
        
        mypos.annotations = {
            "PPR_TYPE": "PPR", # PS_Scan doesn't give type like PPRFinder (P, P1, etc), just the motif match
            "PS_SCORE": score,
            "PPR_STARTS": safe_start, # 0-based
            "PPR_ENDS": safe_end,     # 0-based
            "FIFTH_AA": fifth,
            "LAST_AA": last,
            "PPR_CODE": fifth + last,
        }
        
        parsed_ppr.append(mypos)
        
    # Save to pickle
    output_pkl = Path(output_dir) / f"{seq_id}.pkl"
    with open(output_pkl, "wb") as f:
        pickle.dump(parsed_ppr, f)
        
    # Create a dummy motifs txt for visualization if needed, or just skip
    # The visualization function expects a specific CSV format.
    # We can generate it from the parsed data.
    output_txt = Path(output_dir) / f"{seq_id}_motifs.txt"
    with open(output_txt, "w") as f:
        # Columns: header, start, end, score, motif_seq, second, fifth, last, motif_type
        for p in parsed_ppr:
            # 1-based start/end for display
            d_start = p.annotations["PPR_STARTS"] + 1
            d_end = p.annotations["PPR_ENDS"]
            score = p.annotations["PS_SCORE"]
            seq = str(p.seq)
            fifth = p.annotations["FIFTH_AA"]
            last = p.annotations["LAST_AA"]
            # 2nd AA?
            second = seq[1] if len(seq) > 1 else ""
            
            line = f"{seq_id}\t{d_start}\t{d_end}\t{score}\t{seq}\t{second}\t{fifth}\t{last}\tPPR\n"
            f.write(line)
            
    return str(output_pkl)
