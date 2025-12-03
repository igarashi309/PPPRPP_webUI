import sys
import os
from pathlib import Path

def parse_domtblout(domtblout_path):
    """
    Parses HMMER domtblout file.
    Returns a list of hits: {'name': str, 'start': int, 'end': int, 'score': float, 'evalue': float}
    """
    hits = []
    with open(domtblout_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 22:
                continue
            
            # domtblout format (standard):
            # target_name        accession   tlen query_name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc  description of target
            # 0                  1           2    3                    4           5      6        7      8      9  10  11        12        13     14    15      16  17      18  19      20  21   22...
            
            target_name = parts[0] # Motif name (e.g., P, L1)
            # query_name = parts[3]
            
            # We use domain scores/evalues (cols 11-14 usually, but let's check standard)
            # Col 6: full sequence E-value
            # Col 7: full sequence score
            # Col 11: c-Evalue (conditional)
            # Col 12: i-Evalue (independent)
            # Col 13: domain score
            # Col 14: domain bias
            
            # Coordinates (ali coord): 17 (from), 18 (to) -> 1-based
            # Env coord: 19, 20
            
            try:
                score = float(parts[13])
                evalue = float(parts[12]) # i-Evalue
                start = int(parts[19]) # Using envelope coordinates often better for coverage
                end = int(parts[20])
                
                # Filter weak hits (heuristic similar to original)
                # Original uses E-value < 0.1 or similar?
                # User log showed -E 0.1 in hmmsearch command.
                # We will accept all and filter later if needed.
                
                hits.append({
                    'name': target_name,
                    'start': start,
                    'end': end,
                    'score': score,
                    'evalue': evalue
                })
            except ValueError:
                continue
                
    return hits

def resolve_overlaps(hits):
    """
    Selects the best non-overlapping chain of motifs.
    Strategy: Greedy selection by score.
    """
    # Sort by score descending
    sorted_hits = sorted(hits, key=lambda x: x['score'], reverse=True)
    
    selected_hits = []
    
    for hit in sorted_hits:
        # Check overlap with already selected
        is_overlapping = False
        for selected in selected_hits:
            # Overlap condition: start1 <= end2 and start2 <= end1
            # We allow small overlap? Original PPRFinder allows small overlap?
            # Strict non-overlap for now.
            if hit['start'] <= selected['end'] and selected['start'] <= hit['end']:
                # Calculate overlap amount
                overlap_start = max(hit['start'], selected['start'])
                overlap_end = min(hit['end'], selected['end'])
                overlap_len = overlap_end - overlap_start + 1
                
                # If overlap is significant (> X% of shorter motif?), reject.
                # Let's assume strict for now, or allow very small (e.g. 1-2aa).
                if overlap_len > 5: # Allow 5aa overlap
                    is_overlapping = True
                    break
        
        if not is_overlapping:
            selected_hits.append(hit)
            
    # Sort by position
    selected_hits.sort(key=lambda x: x['start'])
    return selected_hits

def get_sequence_from_fasta(fasta_path):
    """Reads the first sequence from a FASTA file."""
    with open(fasta_path, 'r') as f:
        header = ""
        seq = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if not header:
                    header = line[1:]
                else:
                    break # Only read first sequence
            else:
                seq.append(line)
    return header, "".join(seq)

def write_motifs_txt(hits, output_path, protein_fasta):
    """
    Writes the output in the format expected by parse_pprfinder_output.
    Columns: id, start, end, score, motif_seq, second, fifth_aa, last_aa, motif_type
    """
    header, sequence = get_sequence_from_fasta(protein_fasta)
    
    with open(output_path, 'w') as f:
        # No header in the expected format (pandas read_table header=None)
        
        for hit in hits:
            # Coordinates are 1-based from HMMER
            # Python slicing is 0-based
            # HMMER 'start' is inclusive 1-based. 'end' is inclusive 1-based.
            # So slice is [start-1 : end]
            
            start = hit['start']
            end = hit['end']
            
            # Ensure bounds
            s_idx = max(0, start - 1)
            e_idx = min(len(sequence), end)
            
            motif_seq = sequence[s_idx:e_idx]
            
            # Extract key residues
            # 2nd AA: index 1
            # 5th AA: index 4
            # Last AA: index -1
            
            second = motif_seq[1] if len(motif_seq) > 1 else "-"
            fifth = motif_seq[4] if len(motif_seq) > 4 else "-"
            last = motif_seq[-1] if len(motif_seq) > 0 else "-"
            
            # Columns: id, start, end, score, motif_seq, second, fifth_aa, last_aa, motif_type
            # id seems unused by parser, use protein header or generic
            row = [
                header,
                str(start),
                str(end),
                str(hit['score']),
                motif_seq,
                second,
                fifth,
                last,
                hit['name']
            ]
            
            f.write("\t".join(row) + "\n")

def run(domtblout_path, output_path, protein_fasta):
    hits = parse_domtblout(domtblout_path)
    selected_hits = resolve_overlaps(hits)
    write_motifs_txt(selected_hits, output_path, protein_fasta)
    return output_path

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python ppr_finder.py <domtblout> <output_txt> <protein_fasta>")
        sys.exit(1)
    
    run(sys.argv[1], sys.argv[2], sys.argv[3])
