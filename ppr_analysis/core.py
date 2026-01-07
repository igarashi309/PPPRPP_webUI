"""
Core logic for PPR motif analysis.
Optimized with NumPy.
"""
import os
import numpy as np
from .config import KOBAYASHI_BINDING_SCORES, YAN_BINDING_SCORES, NUCLEOTIDE_INDEX

def dna_to_rna(dna_seq):
    """
    Converts DNA sequence to RNA.
    """
    return dna_seq.replace('T', 'U').replace('t', 'u').upper()

def get_reverse_complement(seq):
    """
    Returns the reverse complement of a DNA/RNA sequence.
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A', 
                  'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'u': 'a',
                  'N': 'N', 'n': 'n'}
    rev_seq = "".join(complement.get(base, base) for base in reversed(seq))
    return rev_seq

def generate_gap_variations(motif_labels, ppr_scores, min_front, min_back, max_gap):
    """
    Generates variations of motif labels with gaps.
    """
    gap_variations = [(motif_labels, ppr_scores, 0, 0)]  # GAP 0, Pos 0

    
    # Optimization: User requested ONLY Gap 0 and Gap N (target).
    # Skip intermediate gaps (1..N-1) to save time.
    target_gaps = [max_gap] if max_gap > 0 else []
    
    for gap_size in target_gaps:
        for front_size in range(min_front, len(motif_labels) - min_back + 1):
            gap_variation = motif_labels[:front_size] + ['--'] * gap_size + motif_labels[front_size:]
            gap_scores = ppr_scores[:front_size] + [0] * gap_size + ppr_scores[front_size:]
            gap_variations.append((gap_variation, gap_scores, gap_size, front_size))
    return gap_variations

def calculate_binding_score_numpy(rna_seq, motif_binding_scores, motif_labels, ppr_scores):
    """
    Calculates binding scores using NumPy for performance.
    """
    # 1. Convert RNA sequence to integer array
    # A=0, C=1, G=2, U=3, Others=-1
    # We can use a lookup table for fast conversion
    # ASCII mapping: A=65, C=67, G=71, U=85
    
    # Create a mapping array (size 256 for ASCII)
    char_to_int = np.full(256, -1, dtype=np.int8)
    char_to_int[ord('A')] = 0
    char_to_int[ord('C')] = 1
    char_to_int[ord('G')] = 2
    char_to_int[ord('U')] = 3
    char_to_int[ord('T')] = 3 # Treat T as U
    
    # Convert sequence string to byte array then to numpy array
    seq_bytes = np.frombuffer(rna_seq.encode('ascii'), dtype=np.uint8)
    seq_indices = char_to_int[seq_bytes]
    
    # 2. Prepare Motif Score Matrix (M x 4)
    # M = number of motifs
    # Rows = motifs, Cols = A, C, G, U scores
    motif_length = len(motif_labels)
    score_matrix = np.zeros((motif_length, 4), dtype=np.float32)
    
    # Base score for unknown/gap is 0.25 (or 0 for gap?)
    # Original logic:
    # if motif == "--": base_score = 0 (from dict)
    # else if motif in dict: base_score = dict[motif][nuc_idx]
    # else: base_score = 0.25
    
    # We also need to multiply by ppr_scores[j] (confidence)
    # score = base_score * ppr_scores[j]
    
    for i, motif in enumerate(motif_labels):
        confidence = ppr_scores[i]
        
        if motif == "--":
            # Gap: scores are 0
            # score_matrix[i, :] = 0 # Already zero
            pass
        elif motif in motif_binding_scores:
            # Valid motif
            base_scores = np.array(motif_binding_scores[motif], dtype=np.float32)
            score_matrix[i, :] = base_scores * confidence
        else:
            # Unknown motif: base 0.25
            score_matrix[i, :] = 0.25 * confidence
            
    # 3. Calculate Scores using Sliding Window
    # We want to compute the sum of scores for each window of size M
    # Input: seq_indices (N)
    # Matrix: score_matrix (M, 4)
    
    # We can use fancy indexing to get scores for each position
    # But we need to handle -1 (invalid chars)
    # Let's replace -1 with 0 (A) temporarily but handle the score separately?
    # Or add a 5th column for invalid chars with score 0?
    # Original logic: if nuc not in NUCLEOTIDE_INDEX: score = 0.25 * confidence (same as unknown motif logic?)
    # Wait, original logic:
    # if motif in motif_binding_scores and nuc in NUCLEOTIDE_INDEX: ...
    # else: score = 0.25 * ppr_scores[j]
    
    # So if nuc is invalid (e.g. 'N'), it gets 0.25 * confidence.
    # Let's add a 5th column for "Other"
    score_matrix_extended = np.zeros((motif_length, 5), dtype=np.float32)
    score_matrix_extended[:, :4] = score_matrix
    
    # Fill 5th column (index 4)
    for i, motif in enumerate(motif_labels):
        confidence = ppr_scores[i]
        if motif == "--":
             score_matrix_extended[i, 4] = 0 # Gap is 0
        else:
             score_matrix_extended[i, 4] = 0.25 * confidence
             
    # Update seq_indices: -1 becomes 4
    seq_indices[seq_indices == -1] = 4
    
    # Now we have seq_indices (N) with values 0-4
    # We want to sum scores over window M
    
    # Construct a matrix of scores for the whole sequence based on motifs?
    # No, the motif matches at relative positions.
    # At position i in sequence, the score is sum(score_matrix[j, seq[i+j]] for j in 0..M-1)
    
    # This is equivalent to:
    # 1. Map entire sequence to scores for EACH motif position? No, that's M x N. Too big.
    # 2. Use sliding window view.
    
    # seq_indices size N
    # Windows: (N - M + 1, M)
    from numpy.lib.stride_tricks import sliding_window_view
    
    if len(seq_indices) < motif_length:
        return []
        
    windows = sliding_window_view(seq_indices, window_shape=motif_length)
    # windows shape: (num_windows, M)
    
    # We need to look up scores.
    # score_matrix_extended shape: (M, 5)
    # We want result shape: (num_windows)
    
    # Advanced indexing:
    # We want result[w] = sum(score_matrix_extended[j, windows[w, j]] for j in 0..M-1)
    
    # We can use `np.take_along_axis` or similar?
    # Or just iterate over M (motif length is small, ~35)?
    # Iterating over M is much faster than iterating over N.
    
    total_scores = np.zeros(len(windows), dtype=np.float32)
    
    # Also need to return detailed info for each site?
    # The original function returns: (i, score, segment, scores, score, motif_binding_score_list)
    # Storing all this for 500k positions is memory intensive and slow to construct.
    # BUT the app likely filters or only displays top hits?
    # Wait, `process_fasta_sequences` returns ALL scores.
    # And `app.py` puts them into a DataFrame.
    # If N=500k, DataFrame has 500k rows. This is why it's slow/memory heavy.
    
    # Optimization:
    # 1. Compute total scores fast.
    # 2. Only construct full details for top hits?
    # But the user might want to see the plot which uses all data?
    # The plot uses `df['Score']`.
    # `visualize_motifs` uses `df`.
    
    # If we must return everything, we can at least optimize the calculation.
    # Constructing the list of tuples is the slow part in Python.
    
    # Let's calculate total scores first.
    for j in range(motif_length):
        # For motif position j, get the nucleotide at window position j
        # windows[:, j] are the nucleotides at pos j for all windows
        col_indices = windows[:, j]
        # Look up scores for motif j
        # score_matrix_extended[j, :] is 1x5
        # We want score_matrix_extended[j, col_indices]
        scores_at_j = score_matrix_extended[j, col_indices]
        total_scores += scores_at_j
        
    # Now we have total scores.
    # To maintain compatibility, we need to return the list of tuples.
    # But generating 500k tuples is slow.
    # Can we return a simplified structure or generator?
    # The `process_fasta_sequences` function consumes this.
    
    # Let's try to optimize the list construction or just return the arrays
    # and let `process_fasta_sequences` handle it?
    # `process_fasta_sequences` iterates and appends to `all_scores`.
    
    # If we change `calculate_binding_score` to return arrays, we need to update `process_fasta_sequences`.
    
    return total_scores, windows, score_matrix_extended

def calculate_binding_score(rna_seq, motif_binding_scores, motif_labels, ppr_scores):
    """
    Wrapper to maintain backward compatibility if needed, 
    but we should update process_fasta_sequences to use the numpy version.
    """
    # This is the old slow version. We replace it or update callers.
    # For this refactor, I will update `process_fasta_sequences` to use `calculate_binding_score_numpy`.
    pass

def process_fasta_sequences(fasta_sequence, motif_binding_scores, motif_labels_variations, ppr_name, predicted_scores, use_ps_score=True):
    """
    Processes FASTA sequence with given variations.
    Optimized to use NumPy.
    """
    all_scores = []
    rna_seq = dna_to_rna(fasta_sequence)

    for motif_labels, ppr_scores, gap_size, gap_pos in motif_labels_variations:
        # Use NumPy calculation
        total_scores, windows, score_matrix = calculate_binding_score_numpy(rna_seq, motif_binding_scores, motif_labels, ppr_scores)
        
        motif_labels_gap = ','.join(motif_labels)
        motif_length = len(motif_labels)
        
        # We need to reconstruct the detailed info for the output
        # But doing it for ALL positions is too slow.
        # Does the app really need 'scores' (list of individual motif scores) and 'motif_binding_score_list' for EVERY position?
        # The CSV output includes them.
        # If we want to support 500kb+, we should probably NOT save lists in CSV for every row.
        # But to keep compatibility, we might have to.
        
        # However, we can optimize the loop.
        # Or better, we can filter? No, user wants "scan".
        
        # Let's iterate and construct. It will still be faster than pure Python calculation + construction.
        # But constructing 500k strings/lists is the bottleneck.
        
        # Maybe we can return a generator?
        # Or just use the total_scores for the plot and only generate details on demand?
        # The `app.py` expects a list of lists/tuples to create a DataFrame.
        
        # Let's try to be as fast as possible in construction.
        
        # Pre-calculate motif_binding_score_list (base scores)?
        # It depends on the sequence.
        
        # Optimization:
        # The return format is:
        # (pos, ppr_binding_score, ppr_name, motif_labels, motif_labels_gap, gap_size, segment, ppr_motif_scores, total_score, motif_binding_score_list, predicted_scores, ppr_scores, gap_pos)
        
        # Most of these are constant or simple.
        # segment, ppr_motif_scores, motif_binding_score_list vary.
        
        # If we just return the data needed for the DataFrame:
        # Position, Score, ...
        
        # Let's iterate efficiently.
        num_sites = len(total_scores)
        
        # Pre-convert constant lists to strings/lists to avoid re-creation if possible?
        # predicted_scores is constant.
        
        # To speed up, we can avoid creating the detailed lists unless necessary.
        # But the CSV writer expects them.
        
        # Let's just implement the loop. It will be faster than before because score calc is done.
        
        # We need to reconstruct individual scores for the output lists.
        # score_matrix (M, 5)
        # windows (N, M) -> indices 0-4
        
        # We can compute all individual scores at once?
        # individual_scores = score_matrix[np.arange(M), windows] -> (N, M)
        # This gives the score for each motif at each position.
        
        # Create indices for broadcasting
        # M = motif_length
        # windows is (N, M)
        # We want to lookup score_matrix[j, windows[i, j]]
        
        # individual_scores = np.take_along_axis(score_matrix, windows, axis=None?) No.
        
        # individual_scores = np.zeros((num_sites, motif_length), dtype=np.float32)
        # for j in range(motif_length):
        #     individual_scores[:, j] = score_matrix[j, windows[:, j]]
            
        # This `individual_scores` corresponds to `scores` in original code.
        # `motif_binding_score_list` is similar but without `ppr_scores` multiplication.
        # We can reconstruct it: individual_scores / ppr_scores[j] ?
        # Or just re-lookup from a base_score_matrix.
        
        # Let's create base_score_matrix
        base_score_matrix = np.zeros((motif_length, 5), dtype=np.float32)
        for i, motif in enumerate(motif_labels):
            if motif == "--":
                 pass
            elif motif in motif_binding_scores:
                base_score_matrix[i, :4] = motif_binding_scores[motif]
            else:
                base_score_matrix[i, :4] = 0.25
            base_score_matrix[i, 4] = 0.25 if motif != "--" else 0
            
        # Compute all individual scores (N x M)
        # This might be memory intensive for 500k x 35 ~ 17.5M floats ~ 70MB. Acceptable.
        
        all_scores_matrix = np.zeros((num_sites, motif_length), dtype=np.float32)
        all_base_scores_matrix = np.zeros((num_sites, motif_length), dtype=np.float32)
        
        for j in range(motif_length):
            col_indices = windows[:, j]
            all_scores_matrix[:, j] = score_matrix[j, col_indices]
            all_base_scores_matrix[:, j] = base_score_matrix[j, col_indices]
            
        # ppr_motif_scores = base_scores * predicted_scores
        # predicted_scores is list of length M?
        # If predicted_scores is empty (default), then ppr_motif_scores is empty?
        # Original code: [b * p for b, p in zip(motif_binding_score_list, predicted_scores)]
        # If predicted_scores is empty, result is empty.
        
        has_predicted = len(predicted_scores) > 0
        if has_predicted:
            # Handle length mismatch if gaps
            # If gap, predicted_scores might need adjustment?
            # Original code just zips, so it truncates.
            # We replicate that.
            limit = min(motif_length, len(predicted_scores))
            pred_array = np.array(predicted_scores[:limit], dtype=np.float32)
            # We only multiply up to limit
            # ppr_motif_scores_matrix = all_base_scores_matrix[:, :limit] * pred_array
            pass
            
        # Construct the columnar data
        # To make it fast, we avoid creating objects for every row if possible.
        
        # Constant columns
        # PPR_Name, Motif_Labels, Motif_Labels_GAP, Gap_Size, Predicted_Scores, PPR_Scores
        # These are same for all rows in this variation.
        
        # We can use numpy repeat/tile if we want arrays, or just lists.
        # Since we are returning to Python, lists might be okay if we don't iterate.
        # But creating a list of 500k identical strings is still memory/time.
        # Pandas handles scalar expansion if we pass a single value? 
        # No, we are concatenating multiple variations (gap sizes).
        # So we must expand them to match num_sites.
        
        # Let's create lists efficiently.
        # [value] * num_sites is fast.
        
        ppr_name_col = [ppr_name] * num_sites
        motif_labels_col = [','.join(motif_labels)] * num_sites
        motif_labels_gap_col = [motif_labels_gap] * num_sites
        gap_size_col = [gap_size] * num_sites
        gap_pos_col = [gap_pos] * num_sites
        predicted_scores_col = [predicted_scores] * num_sites
        ppr_scores_col = [ppr_scores] * num_sites
        
        # Variable columns
        # Position: 0..num_sites-1
        pos_col = list(range(num_sites))
        
        # PPR_Binding_Score (Total Score)
        # total_scores is numpy array.
        ppr_binding_score_col = total_scores.tolist()
        total_score_col = ppr_binding_score_col # Same
        
        if not use_ps_score:
             # Recompute total if needed (sum of base scores)
             # all_base_scores_matrix sum along axis 1
             ppr_binding_score_col = np.sum(all_base_scores_matrix, axis=1).tolist()
             total_score_col = ppr_binding_score_col
             
        # Segment
        # This is the slow part: creating 500k substrings.
        # But we can't avoid it if the user needs to see the sequence.
        # List comprehension is fastest way in Python.
        segment_col = [rna_seq[i:i+motif_length] for i in range(num_sites)]
        
        # PPR_Motif_Scores (List of floats)
        # This is also slow: creating 500k lists.
        # all_scores_matrix is (N, M).
        # We want list of lists.
        # all_scores_matrix.tolist() is optimized in C.
        # But we need to handle predicted_scores multiplication if present.
        
        if has_predicted:
            limit = min(motif_length, len(predicted_scores))
            # Multiply relevant columns
            # We need to use base scores for this multiplication
            weighted_scores = all_base_scores_matrix[:, :limit] * np.array(predicted_scores[:limit], dtype=np.float32)
            # We need to append the rest of the scores (if any) or just return weighted?
            # Original logic: zip(motif_binding_score_list, predicted_scores)
            # So result length is min(len(motifs), len(predicted)).
            ppr_motif_scores_col = weighted_scores.tolist()
        else:
            ppr_motif_scores_col = [[]] * num_sites
            
        # PPR_Binding_Scores_List (Base scores)
        motif_binding_score_list_col = all_base_scores_matrix.tolist()
        
        # Append to main lists
        # We are aggregating across variations (gap sizes).
        # So we extend the lists.
        
        if not all_scores:
            # Initialize dict
            all_scores = {
                "Position_Relative": [],
                "PPR_Binding_Score": [],
                "PPR_Name": [],
                "Motif_Labels": [],
                "Motif_Labels_GAP": [],
                "Gap_Size": [],
                "Segment": [],
                "PPR_Motif_Scores": [],
                "Total_Score": [],
                "PPR_Binding_Scores_List": [],
                "Predicted PPR Motif Scores List": [],
                "Motif_Confidence_Scores": [],
                "Gap_Index": [] # Add new column
            }
            
        all_scores["Position_Relative"].extend(pos_col)
        all_scores["PPR_Binding_Score"].extend(ppr_binding_score_col)
        all_scores["PPR_Name"].extend(ppr_name_col)
        all_scores["Motif_Labels"].extend(motif_labels_col)
        all_scores["Motif_Labels_GAP"].extend(motif_labels_gap_col)
        all_scores["Gap_Size"].extend(gap_size_col)
        all_scores["Segment"].extend(segment_col)
        all_scores["PPR_Motif_Scores"].extend(ppr_motif_scores_col)
        all_scores["Total_Score"].extend(total_score_col)
        all_scores["PPR_Binding_Scores_List"].extend(motif_binding_score_list_col)
        all_scores["Predicted PPR Motif Scores List"].extend(predicted_scores_col)
        all_scores["Motif_Confidence_Scores"].extend(ppr_scores_col)
        all_scores["Gap_Index"].extend(gap_pos_col)

    return all_scores

def analyze_ppr_data(motif_data, fasta_sequence, min_front, min_back, max_gap, predicted_scores, ppr_name):
    """
    Analyzes a single PPR data object (extracted from pickle).
    
    Args:
        motif_data (dict): {'motif_labels': [], 'ppr_scores': []}
        fasta_sequence (str): The sequence.
        min_front, min_back, max_gap (int): Gap parameters.
        predicted_scores (list): List of predicted scores.
        ppr_name (str): Name of the PPR.
        
    Returns:
        dict: Scores for all analysis options.
    """
    motif_labels = motif_data['motif_labels']
    ppr_scores = motif_data['ppr_scores']
    
    gap0_variations = [(motif_labels, ppr_scores, 0, 0)]
    gap_variations = generate_gap_variations(motif_labels, ppr_scores, min_front, min_back, max_gap)
    
    scores = {
        'kobayashi_gap0_ps': process_fasta_sequences(fasta_sequence, KOBAYASHI_BINDING_SCORES, gap0_variations, ppr_name, predicted_scores, use_ps_score=True),
        'kobayashi_gap_ps': process_fasta_sequences(fasta_sequence, KOBAYASHI_BINDING_SCORES, gap_variations, ppr_name, predicted_scores, use_ps_score=True),
        'kobayashi_gap0_no_ps': process_fasta_sequences(fasta_sequence, KOBAYASHI_BINDING_SCORES, gap0_variations, ppr_name, predicted_scores, use_ps_score=False),
        'kobayashi_gap_no_ps': process_fasta_sequences(fasta_sequence, KOBAYASHI_BINDING_SCORES, gap_variations, ppr_name, predicted_scores, use_ps_score=False),
        'yan_gap0_ps': process_fasta_sequences(fasta_sequence, YAN_BINDING_SCORES, gap0_variations, ppr_name, predicted_scores, use_ps_score=True),
        'yan_gap_ps': process_fasta_sequences(fasta_sequence, YAN_BINDING_SCORES, gap_variations, ppr_name, predicted_scores, use_ps_score=True),
        'yan_gap0_no_ps': process_fasta_sequences(fasta_sequence, YAN_BINDING_SCORES, gap0_variations, ppr_name, predicted_scores, use_ps_score=False),
        'yan_gap_no_ps': process_fasta_sequences(fasta_sequence, YAN_BINDING_SCORES, gap_variations, ppr_name, predicted_scores, use_ps_score=False)
    }
    return scores

