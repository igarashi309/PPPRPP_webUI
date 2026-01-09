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
    # Skip intermediate gaps (1..N-1) to save time, as confirmed by User.
    target_gaps = [int(max_gap)] if max_gap > 0 else []
    
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
    all_scores = {
        "Position_Relative": [],
        "PPR_Binding_Score": [],
        "PPR_Name": [],
        "Motif_Labels": [],
        "Motif_Labels_GAP": [],
        "Gap_Size": [],
        "Gap_Index": [],
        "Segment": [],
        "PPR_Motif_Scores": [],
        "Total_Score": [],
        "PPR_Binding_Scores_List": [],
        "Predicted PPR Motif Scores List": [],
        "Motif_Confidence_Scores": []
    }
    
    rna_seq = dna_to_rna(fasta_sequence)
    seq_len = len(rna_seq)
    
    # Pre-calculate Total Scores for ALL variations
    # Result: list of (total_scores_array, variation_obj)
    variation_results = []
    
    for v_idx, variation in enumerate(motif_labels_variations):
        motif_labels, ppr_scores_original, gap_size, gap_pos = variation
        
        # Adjust ppr_scores for calculation based on use_ps_score
        ppr_scores_for_calc = ppr_scores_original
        if not use_ps_score:
            # If not using PS scores, treat all confidence scores as 1.0 for the binding score calculation
            ppr_scores_for_calc = [1.0] * len(motif_labels)

        # Use NumPy calculation
        # We only need total_scores first to find the best variation
        total_scores, windows, score_matrix = calculate_binding_score_numpy(rna_seq, motif_binding_scores, motif_labels, ppr_scores_for_calc)
        
        # total_scores length might vary if window size varies (Gap 0 vs Gap 5)
        # But for 'scan', results are aligned to start position?
        # sliding_window_view returns N - M + 1 windows.
        # If M varies (Gap 0 len 30 vs Gap 5 len 35), result length differs by 5.
        # This makes alignment tricky.
        # We need to align by "Position_Relative" which is index in sequence.
        # A hit at pos 10 means the motif starts at 10.
        # If L=30, it covers 10-39. If L=35, it covers 10-44.
        # But both are "at pos 10".
        # So we should pad/truncate to match sequence length?
        # Or just handle the common overlap.
        # The number of valid windows is `seq_len - motif_length + 1`.
        # Max valid index is `seq_len - motif_length`.
        # So Gap 0 (shorter) has MORE valid positions than Gap 5 (longer).
        # We can only compare positions valid for BOTH (or ALL).
        # i.e., up to `seq_len - max_len`.
        
        variation_results.append({
            "total_scores": total_scores,
            "windows": windows,
            "score_matrix": score_matrix,
            "variation": variation, # Store original variation tuple
            "motif_len": len(motif_labels)
        })

    if not variation_results:
        return all_scores

    # Find common valid length (limited by longest motif)
    max_motif_len = max(v["motif_len"] for v in variation_results)
    num_valid_sites = seq_len - max_motif_len + 1
    
    if num_valid_sites <= 0:
        return all_scores

    # Stack total scores for comparison
    # We truncate all score arrays to num_valid_sites
    # (Shorter motifs have more scores, we ignore the tail where longer motifs don't fit)
    
    stacked_scores = np.zeros((len(variation_results), num_valid_sites), dtype=np.float32)
    
    for i, res in enumerate(variation_results):
        scores = res["total_scores"]
        # If scores is shorter than num_valid_sites (shouldn't happen if logic correct), pad with -inf?
        # Actual logic: valid_sites = len - len(motif) + 1.
        # Longest motif has Shortest scores array.
        # So `num_valid_sites` IS the length of the shortest scores array.
        # All other arrays are LONGER. We truncate them.
        if len(scores) < num_valid_sites:
            # This should ideally not happen if num_valid_sites is calculated correctly
            # based on max_motif_len. If it does, it means a motif is too long for the sequence.
            # Fill with a very low score so it's not chosen.
            stacked_scores[i, :len(scores)] = scores
            stacked_scores[i, len(scores):] = -np.inf
        else:
            stacked_scores[i, :] = scores[:num_valid_sites]
        
    # Find best variation for each position
    best_indices = np.argmax(stacked_scores, axis=0) # Shape (num_valid_sites,)
    
    # Now construct the result columns
    # We create lists of length `num_valid_sites`
    
    # Initialize output lists
    out_pos = list(range(num_valid_sites))
    # We need to fill other lists based on best_indices
    
    # To optimize, we loop over variations and fill the positions where that variation won
    
    # Pre-allocate lists with None or placeholders
    # (Using None allows us to spot bugs if holes remain)
    out_ppr_name = [None] * num_valid_sites
    out_motif_labels = [None] * num_valid_sites
    out_motif_labels_gap = [None] * num_valid_sites
    out_gap_size = [None] * num_valid_sites
    out_gap_index = [None] * num_valid_sites
    out_segment = [None] * num_valid_sites
    out_ppr_motif_scores = [None] * num_valid_sites
    out_total_score = [None] * num_valid_sites
    out_binding_scores_list = [None] * num_valid_sites
    out_pred_scores_list = [None] * num_valid_sites
    out_conf_scores = [None] * num_valid_sites
    
    for v_idx, res in enumerate(variation_results):
        # Identify positions where this variation is best
        mask = (best_indices == v_idx)
        if not np.any(mask):
            continue
            
        indices = np.where(mask)[0]
        
        # Extract data for this variation
        vars_tuple = res["variation"]
        motif_labels, ppr_scores_original, gap_size, gap_pos = vars_tuple
        windows = res["windows"]
        score_matrix = res["score_matrix"] # This score_matrix was built using ppr_scores_for_calc
        
        # Constant strings for this variation
        m_labels_str = ','.join(motif_labels)
        m_labels_gap_str = ','.join(motif_labels) # Same as above but variable name preserved
        
        # Fill simple scalar/string columns
        for idx in indices:
            out_ppr_name[idx] = ppr_name
            out_motif_labels[idx] = m_labels_str
            out_motif_labels_gap[idx] = m_labels_gap_str # Logic same as original
            out_gap_size[idx] = gap_size
            out_gap_index[idx] = gap_pos
            out_pred_scores_list[idx] = predicted_scores
            out_conf_scores[idx] = ppr_scores_original # Use original ppr_scores for output
            out_total_score[idx] = float(stacked_scores[v_idx, idx])
            
        # For complex columns (Segment, detailed scores), we calculate only for needed indices
        
        # Base Score Matrix (without ppr_scores multiplication)
        motif_length = len(motif_labels)
        base_score_matrix = np.zeros((motif_length, 5), dtype=np.float32)
        for i, motif in enumerate(motif_labels):
            if motif == "--":
                 pass
            elif motif in motif_binding_scores:
                base_score_matrix[i, :4] = motif_binding_scores[motif]
            else:
                base_score_matrix[i, :4] = 0.25
            base_score_matrix[i, 4] = 0.25 if motif != "--" else 0
            
        # Predicted Scores logic
        has_predicted = len(predicted_scores) > 0
        pred_array = None
        if has_predicted:
             limit = min(motif_length, len(predicted_scores))
             pred_array = np.array(predicted_scores[:limit], dtype=np.float32)
        
        # Calculate details for winning indices
        for idx in indices:
            # Segment
            out_segment[idx] = rna_seq[idx : idx + motif_length]
            
            # Scores
            # windows[idx] is the column indices for this position
            # Note: windows array was calculated for full length, idx is valid.
            col_indices = windows[idx] 
            
            seq_range = np.arange(motif_length)
            
            # scores_vec uses ppr_scores_for_calc (which might be all 1s if use_ps_score=False)
            scores_vec = score_matrix[seq_range, col_indices]
            base_scores_vec = base_score_matrix[seq_range, col_indices]
            
            out_binding_scores_list[idx] = base_scores_vec.tolist()
            
            if has_predicted:
                limit = min(motif_length, len(predicted_scores))
                weighted = base_scores_vec[:limit] * pred_array
                out_ppr_motif_scores[idx] = weighted.tolist()
            else:
                out_ppr_motif_scores[idx] = []
                
    # Final assembly
    all_scores["Position_Relative"] = out_pos
    all_scores["PPR_Name"] = out_ppr_name
    all_scores["Motif_Labels"] = out_motif_labels
    all_scores["Motif_Labels_GAP"] = out_motif_labels_gap
    all_scores["Gap_Size"] = out_gap_size
    all_scores["Gap_Index"] = out_gap_index
    all_scores["Segment"] = out_segment
    all_scores["PPR_Motif_Scores"] = out_ppr_motif_scores
    all_scores["Total_Score"] = out_total_score
    all_scores["PPR_Binding_Scores_List"] = out_binding_scores_list
    all_scores["Predicted PPR Motif Scores List"] = out_pred_scores_list
    all_scores["Motif_Confidence_Scores"] = out_conf_scores
    
    # Handle Total Score override if no_ps
    # If use_ps_score is False, PPR_Binding_Score should be sum of base scores, not weighted scores.
    if not use_ps_score:
        # Re-sum base scores for the PPR_Binding_Score column
        ppr_binding_score_col = []
        for i in range(num_valid_sites):
             if out_binding_scores_list[i]:
                 ppr_binding_score_col.append(sum(out_binding_scores_list[i]))
             else:
                 ppr_binding_score_col.append(0.0) # Ensure float
        all_scores["PPR_Binding_Score"] = ppr_binding_score_col
    else:
        # If use_ps_score is True, PPR_Binding_Score is the same as Total_Score (which is weighted)
        all_scores["PPR_Binding_Score"] = out_total_score

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

def analyze_single_site(motif_data, rna_seq, site_pos, min_front, min_back, max_gap, predicted_scores, use_ps_score, method):
    """
    Analyzes a single genomic site for all gap variations.
    Used for "Detail Inspector" to show alternative alignments.
    """
    import pandas as pd
    import numpy as np

    # Determine scoring table
    if method == 'kobayashi':
         motif_binding_scores = KOBAYASHI_BINDING_SCORES
    else:
         motif_binding_scores = YAN_BINDING_SCORES
         
    # Generate variations
    variations = generate_gap_variations(motif_data['motif_labels'], motif_data['ppr_scores'], min_front, min_back, max_gap)
    
    results = []
    
    # Adjust ppr_scores based on use_ps_score
    ppr_scores_for_calc = motif_data['ppr_scores']
    if not use_ps_score:
        ppr_scores_for_calc = [1.0] * len(motif_data['motif_labels'])
        
    # Re-generate variations with correct calc scores
    calc_variations = generate_gap_variations(motif_data['motif_labels'], ppr_scores_for_calc, min_front, min_back, max_gap)
    
    seq_len = len(rna_seq)
    
    # Iterate both lists (they are aligned)
    for (v_orig, v_calc) in zip(variations, calc_variations):
        # v_orig: (labels, original_scores, gap_size, front_size)
        # v_calc: (labels, calc_scores, gap_size, front_size)
        
        m_labels = v_calc[0]
        m_scores_calc = v_calc[1]
        gap_size = v_calc[2]
        front_size = v_calc[3] # Index of gap start
        
        motif_len = len(m_labels)
        
        # Check bounds
        if site_pos < 0 or site_pos + motif_len > seq_len:
            continue
            
        # Extract sub-sequence
        sub_seq = rna_seq[site_pos : site_pos + motif_len]
        
        # Calculate score for this single window
        # calculate_binding_score_numpy takes seq.
        total_scores, windows, score_matrix = calculate_binding_score_numpy(sub_seq, motif_binding_scores, m_labels, m_scores_calc)
        
        if len(total_scores) == 0:
            continue
            
        score = float(total_scores[0])
        
        # Base scores for detail display
        base_score_matrix = np.zeros((motif_len, 5), dtype=np.float32)
        for i, motif in enumerate(m_labels):
            if motif == "--": pass
            elif motif in motif_binding_scores: base_score_matrix[i, :4] = motif_binding_scores[motif]
            else: base_score_matrix[i, :4] = 0.25
            base_score_matrix[i, 4] = 0.25 if motif != "--" else 0
            
        # Get indices for this window (windows[0])
        col_indices = windows[0]
        seq_range = np.arange(motif_len)
        
        base_scores_vec = base_score_matrix[seq_range, col_indices]
        
        # Original confidence scores (for display)
        conf_scores_gap = v_orig[1]
        
        # Predicted scores (placeholder/basic logic)
        pred_scores_list = predicted_scores if predicted_scores else []
        ppr_motif_scores = [] 
        if pred_scores_list:
             limit = min(motif_len, len(pred_scores_list))
             weighted = base_scores_vec[:limit] * np.array(pred_scores_list[:limit], dtype=np.float32)
             ppr_motif_scores = weighted.tolist()
             
        # Gap Index in dataframe terminology (position of gap in motif)
        # Is Gap_Index used for anything? Logic says "Gap_Index" is front_size.
        
        row = {
            "Gap_Index": front_size,
            "Gap_Size": gap_size,
            "Total_Score": score, # Weighted if use_ps_score
            "PPR_Binding_Score": score if use_ps_score else float(np.sum(base_scores_vec)),
            "Motif_Labels_GAP": ','.join(m_labels),
            "Segment": sub_seq,
            "PPR_Binding_Scores_List": base_scores_vec.tolist(),
            "Motif_Confidence_Scores": conf_scores_gap,
            "PPR_Motif_Scores": ppr_motif_scores,
            "Predicted PPR Motif Scores List": pred_scores_list
        }
        results.append(row)
        
    return pd.DataFrame(results)

