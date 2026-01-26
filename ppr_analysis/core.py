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

def generate_gap_variations(motif_labels, ppr_scores, min_front, min_back, max_gap, all_gaps=False, motif_types=None):
    """
    Generates variations of motif labels with gaps.
    If all_gaps is True, generates all gap sizes from 1 to max_gap.
    Otherwise (default), generates only gap 0 and max_gap (optimization).
    Includes motif_types variation if provided.
    """
    # Base variation
    # Result tuple: (labels, scores, gap_size, front_size, types_variation)
    
    # If motif_types not provided, fill with empty strings or copy labels? 
    # Better to keep None if not used, but tuple size must be consistent.
    # Let's assume list of types valid if provided.
    
    if motif_types is None:
        motif_types = motif_labels # Fallback/Duplicate
        
    base_var = (motif_labels, ppr_scores, 0, 0, motif_types)
    gap_variations = [base_var]  # GAP 0, Pos 0

    target_gaps = []
    if max_gap > 0:
        if all_gaps:
             target_gaps = list(range(1, int(max_gap) + 1))
        else:
             target_gaps = [int(max_gap)]
    
    for gap_size in target_gaps:
        # Sliding Gap Insertion:
        # Insert ONE block of gap_size at variable positions (front_size).
        
        # Valid positions for gap insertion: 
        # min_front to len(motif_labels) - min_back
        # e.g., if min_front=1, gap can happen after 1st motif.
        
        for front_size in range(min_front, len(motif_labels) - min_back + 1):
             gap_variation = motif_labels[:front_size] + ['--'] * gap_size + motif_labels[front_size:]
             gap_scores = ppr_scores[:front_size] + [0] * gap_size + ppr_scores[front_size:]
             
             # Also gap for types
             gap_types = motif_types[:front_size] + ['--'] * gap_size + motif_types[front_size:]
             
             gap_variations.append((gap_variation, gap_scores, gap_size, front_size, gap_types))
            
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
        return np.array([]), np.array([]), np.array([])
        
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

def process_fasta_sequences(fasta_sequence, motif_binding_scores, motif_labels_variations, ppr_name, predicted_scores, use_ps_score=True, detailed_output=True):
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
        "Motif_Types_GAP": [], # New Column
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
    
    # Optimization: If Lite Mode, we can optimize the loop?
    # No, we need total scores for filter. 
    
    for v_idx, variation in enumerate(motif_labels_variations):
        motif_labels, ppr_scores_original, gap_size, gap_pos, gap_types = variation
        
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

    # Aggregation: Return ALL variations instead of just the best one.
    # We append results from each variation sequentially.
    
    # Initialize output lists (accumulators)
    out_pos = []
    out_ppr_name = []
    out_motif_labels = []
    out_motif_labels_gap = []
    out_motif_types_gap = []
    out_gap_size = []
    out_gap_index = []
    out_segment = []
    out_ppr_motif_scores = []
    out_total_score = []
    out_binding_scores_list = []
    out_pred_scores_list = []
    out_conf_scores = []
    
    # Re-verify scoring method for correct binding score output
    # If no_ps, total_score logic differs? 
    # Logic: PPR_Binding_Score column is used by UI.
    # If use_ps_score=True, Matrix includes confidence. total_scores is Weighted.
    # If use_ps_score=False, Matrix uses 1.0 confidence. total_scores is Unweighted.
    # So `res["total_scores"]` is strictly what we want for "PPR_Binding_Score" column in both cases?
    # Wait, existing code had a check at line 424 to Re-sum?
    # Ah, if use_ps_score=False, `ppr_scores_for_calc` was set to all 1s.
    # So `total_scores` IS the unweighted sum. 
    # So we don't need to re-sum manually if we trust the loop logic. 
    # But let's check: "PPR_Binding_Score: score if use_ps_score else float(np.sum(base_scores_vec))" in analyze_single_site.
    # In `process_fasta_sequences`:
    # Line 245: if not use_ps_score: ppr_scores_for_calc = [1.0]...
    # Line 249: calculate...
    # So `total_scores` ALREADY follows the rule.
    # Why did line 424 exists? "Handle Total Score override if no_ps".
    # It recalculated `ppr_binding_score_col` from `out_binding_scores_list`.
    # `out_binding_scores_list` came from `base_score_matrix` (unweighted).
    # `score_matrix` (used for total_scores) was `base * ppr_scores_for_calc`.
    # If no_ps, `ppr_scores_for_calc` is 1s. So `score_matrix` == `base_score_matrix`.
    # So `total_scores` == `sum(binding_scores_list)`.
    # So the re-calculation was redundant (or verification). 
    # I will rely on `total_scores`.
    
    for v_idx, res in enumerate(variation_results):
        scores = res["total_scores"]
        num_sites = len(scores)
        
        if num_sites == 0:
            continue
            
        vars_tuple = res["variation"]
        motif_labels, ppr_scores_original, gap_size, gap_pos, gap_types = vars_tuple
        windows = res["windows"]
        score_matrix = res["score_matrix"]
        
        # Constant strings
        m_labels_str = ','.join(motif_labels)
        m_types_str = ','.join(gap_types) if gap_types else ""
        
        # Extend constant columns (ALWAYS NEEDED)
        out_pos.extend(range(num_sites))
        out_ppr_name.extend([ppr_name] * num_sites)
        out_gap_size.extend([gap_size] * num_sites)
        out_gap_index.extend([gap_pos] * num_sites)
        out_total_score.extend(scores.tolist())
        
        # Complex columns Loop
        if not detailed_output:
            # Lite Mode: Drop ALL redundant columns that can be inferred from context
            out_segment.extend([""] * num_sites)
            
            # Empty lists / placeholders for heavy columns
            out_binding_scores_list.extend([[]] * num_sites) 
            out_ppr_motif_scores.extend([[]] * num_sites)
            out_pred_scores_list.extend([[]] * num_sites)
            
            # Also drop constant columns repeated millions of times
            out_motif_labels.extend([""] * num_sites) 
            out_motif_labels_gap.extend([""] * num_sites)
            out_motif_types_gap.extend([""] * num_sites)
            out_conf_scores.extend([[]] * num_sites)
            continue
            
        # Detailed Mode: Extend with actual large data
        out_motif_labels.extend([m_labels_str] * num_sites)
        out_motif_labels_gap.extend([m_labels_str] * num_sites)
        out_motif_types_gap.extend([m_types_str] * num_sites)
        out_pred_scores_list.extend([predicted_scores] * num_sites)
        out_conf_scores.extend([ppr_scores_original] * num_sites)
        
        motif_length = len(motif_labels)
        
        # 1. Segments
        # Vectorized segment extraction is hard in strings. List comp is fast enough.
        # rna_seq is string.
        # segments = [rna_seq[i : i + motif_length] for i in range(num_sites)]
        out_segment.extend([rna_seq[i : i + motif_length] for i in range(num_sites)])
        
        # 2. Detailed Scores
        # score_matrix (M, 5). windows (N, M).
        # We want to extract [M] scores for each N.
        # Using fancy indexing for bulk extraction:
        # result shape (N, M).
        # We need to iterate M columns of windows?
        # Re-using logic from calculate_binding_score_numpy or simply iterate?
        # To get list of lists [[s1, s2...], [s1, s2...]], iteration is unavoidable for storage.
        
        # Reconstruction of base_score_matrix (Unweighted) for display
        base_score_matrix = np.zeros((motif_length, 5), dtype=np.float32)
        for i, motif in enumerate(motif_labels):
            if motif == "--": pass
            elif motif in motif_binding_scores: base_score_matrix[i, :4] = motif_binding_scores[motif]
            else: base_score_matrix[i, :4] = 0.25
            base_score_matrix[i, 4] = 0.25 if motif != "--" else 0
            
        # Get all scores at once?
        # windows is (N, M). 
        # We want to map each element of windows to a score from base_score_matrix.
        # base_score_matrix has shape (M, 5). 
        # Element at windows[i, j] (nucleotide index) corresponds to motif pos j.
        # So we want score = base_score_matrix[j, windows[i, j]].
        
        # This can be done via `np.take_along_axis`? No, broadcasting.
        # scores_block = base_score_matrix[np.arange(motif_length), windows] -> No.
        
        # Correct bulk indexing:
        # We iterate cols j=0..M-1.
        # For col j, we take base_score_matrix[j, windows[:, j]].
        # This gives a column of scores (length N).
        # We stack them to get (N, M).
        
        list_of_score_cols = []
        for j in range(motif_length):
            col_indices = windows[:, j]
            scores_col = base_score_matrix[j, col_indices]
            list_of_score_cols.append(scores_col)
            
        # Stack columns to get (N, M)
        # scores_n_m = np.column_stack(list_of_score_cols) # Shape (N, M)
        scores_n_m = np.stack(list_of_score_cols, axis=1)
        
        # Convert to list of lists
        out_binding_scores_list.extend(scores_n_m.tolist())
        
        # Pred scores / Confidence used for Weighting
        # We prefer `ppr_scores_original` from the variation if meaningful, 
        # or `predicted_scores` passed in (often None/Empty in unified analysis if handled by object).
        # Actually `ppr_scores_original` IS the source of truth for this variation.
        
        vector_for_weighting = np.array(ppr_scores_original, dtype=np.float32)
        
        # If predicted_scores is passed explicitly, it might override? 
        # But `ppr_scores_original` comes FROM `motif_data` which comes FROM `ppr_finder`.
        # So `ppr_scores_original` is correct.
        
        # Check alignment
        if len(vector_for_weighting) == motif_length:
             # Broadcast multiply: (N, M) * (M,)
             # We need to reshape vector to (1, M) to broadcast across N rows?
             # scores_n_m is (N, M). vector is (M,).
             # NumPy broadcasts last dimension automatically.
             
             weighted_nm = scores_n_m * vector_for_weighting
             out_ppr_motif_scores.extend(weighted_nm.tolist())
        else:
             # Mismatch length (rare/shouldn't happen for valid variation)
             # Fallback to unweighted or empty?
             # Try unweighted but log warning?
             out_ppr_motif_scores.extend(scores_n_m.tolist()) 
             
        # Maintain predicted_scores column just in case (though likely empty list)
        out_pred_scores_list.extend([predicted_scores] * num_sites)
    
    # Assembly
    all_scores["Position_Relative"] = out_pos
    all_scores["PPR_Name"] = out_ppr_name
    all_scores["Motif_Labels"] = out_motif_labels
    all_scores["Motif_Labels_GAP"] = out_motif_labels_gap
    all_scores["Motif_Types_GAP"] = out_motif_types_gap
    all_scores["Gap_Size"] = out_gap_size
    all_scores["Gap_Index"] = out_gap_index
    all_scores["Segment"] = out_segment
    all_scores["PPR_Motif_Scores"] = out_ppr_motif_scores
    all_scores["Total_Score"] = out_total_score
    all_scores["PPR_Binding_Scores_List"] = out_binding_scores_list
    all_scores["Predicted PPR Motif Scores List"] = out_pred_scores_list
    all_scores["Motif_Confidence_Scores"] = out_conf_scores
    all_scores["PPR_Binding_Score"] = out_total_score

    return all_scores

def analyze_ppr_data(motif_data, fasta_sequence, min_front, min_back, max_gap, predicted_scores, ppr_name, detailed_output=True, target_analysis=None):
    """
    Analyzes PPR data. 
    If target_analysis is provided (dict with keys 'method', 'use_ps', 'gap_mode'), 
    computes ONLY that specific analysis. Otherwise computes all 8 variations (legacy behavior).
    
    target_analysis example: {'method': 'yan', 'use_ps': False, 'gap_mode': 'varied'} 
    gap_mode: 'gap0' or 'varied'
    """
    motif_labels = motif_data['motif_labels']
    ppr_scores = motif_data['ppr_scores']
    motif_types = motif_data.get('motif_types', []) # Get types if available
    if not motif_types:
        motif_types = motif_labels # Fallback
    
    scores = {}
    
    # Helper to determine which to run
    def should_run(method, use_ps, mode):
        if not target_analysis: return True # Run all
        if target_analysis['method'] != method: return False
        if target_analysis['use_ps'] != use_ps: return False
        if target_analysis['gap_mode'] != mode: return False
        return True
        
    gap0_vars = None
    gapN_vars = None
    
    # Kobayashi
    if should_run('kobayashi', True, 'gap0') or should_run('kobayashi', False, 'gap0'):
        gap0_vars = [(motif_labels, ppr_scores, 0, 0, motif_types)]
    if should_run('kobayashi', True, 'varied') or should_run('kobayashi', False, 'varied'):
        gapN_vars = generate_gap_variations(motif_labels, ppr_scores, min_front, min_back, max_gap, motif_types=motif_types)

    # Yan
    if should_run('yan', True, 'gap0') or should_run('yan', False, 'gap0'):
        gap0_vars = [(motif_labels, ppr_scores, 0, 0, motif_types)]
    if should_run('yan', True, 'varied') or should_run('yan', False, 'varied'):
        if gapN_vars is None: gapN_vars = generate_gap_variations(motif_labels, ppr_scores, min_front, min_back, max_gap, motif_types=motif_types)

    # Execute
    if should_run('kobayashi', True, 'gap0'):
        scores['kobayashi_gap0_ps'] = process_fasta_sequences(fasta_sequence, KOBAYASHI_BINDING_SCORES, gap0_vars, ppr_name, predicted_scores, True, detailed_output)
    if should_run('kobayashi', True, 'varied'):
        scores['kobayashi_gap_ps'] = process_fasta_sequences(fasta_sequence, KOBAYASHI_BINDING_SCORES, gapN_vars, ppr_name, predicted_scores, True, detailed_output)
    if should_run('kobayashi', False, 'gap0'):
        scores['kobayashi_gap0_no_ps'] = process_fasta_sequences(fasta_sequence, KOBAYASHI_BINDING_SCORES, gap0_vars, ppr_name, predicted_scores, False, detailed_output)
    if should_run('kobayashi', False, 'varied'):
        scores['kobayashi_gap_no_ps'] = process_fasta_sequences(fasta_sequence, KOBAYASHI_BINDING_SCORES, gapN_vars, ppr_name, predicted_scores, False, detailed_output)
        
    if should_run('yan', True, 'gap0'):
        scores['yan_gap0_ps'] = process_fasta_sequences(fasta_sequence, YAN_BINDING_SCORES, gap0_vars, ppr_name, predicted_scores, True, detailed_output)
    if should_run('yan', True, 'varied'):
        scores['yan_gap_ps'] = process_fasta_sequences(fasta_sequence, YAN_BINDING_SCORES, gapN_vars, ppr_name, predicted_scores, True, detailed_output)
    if should_run('yan', False, 'gap0'):
        scores['yan_gap0_no_ps'] = process_fasta_sequences(fasta_sequence, YAN_BINDING_SCORES, gap0_vars, ppr_name, predicted_scores, False, detailed_output)
    if should_run('yan', False, 'varied'):
        scores['yan_gap_no_ps'] = process_fasta_sequences(fasta_sequence, YAN_BINDING_SCORES, gapN_vars, ppr_name, predicted_scores, False, detailed_output)

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
    # Pass motif_types if available
    motif_types = motif_data.get('motif_types', [])
    if not motif_types:
        motif_types = motif_data['motif_labels']
        
    variations = generate_gap_variations(motif_data['motif_labels'], motif_data['ppr_scores'], min_front, min_back, max_gap, all_gaps=False, motif_types=motif_types)
    
    results = []
    
    # Adjust ppr_scores based on use_ps_score
    ppr_scores_for_calc = motif_data['ppr_scores']
    if not use_ps_score:
        ppr_scores_for_calc = [1.0] * len(motif_data['motif_labels'])
        
    # Re-generate variations with correct calc scores
    calc_variations = generate_gap_variations(motif_data['motif_labels'], ppr_scores_for_calc, min_front, min_back, max_gap, all_gaps=False, motif_types=motif_types)
    
    seq_len = len(rna_seq)
    
    # Iterate both lists (they are aligned)
    for (v_orig, v_calc) in zip(variations, calc_variations):
        # v_orig: (labels, original_scores, gap_size, front_size, types_variation)
        # v_calc: (labels, calc_scores, gap_size, front_size, types_variation)
        
        m_labels = v_calc[0]
        m_scores_calc = v_calc[1]
        gap_size = v_calc[2]
        front_size = v_calc[3] # Index of gap start
        m_types = v_calc[4] # Types
        
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
        
        # Apply weighting from m_scores_calc (contains confidence or 1.0 depending on use_ps_score)
        # m_scores_calc corresponds to motifs in m_labels.
        # We need to map it to the sequence length? No, m_scores_calc is length M.
        # base_scores_vec is length M.
        
        weighted_scores_vec = base_scores_vec * np.array(m_scores_calc, dtype=np.float32)
        
        # Original confidence scores (for display)
        conf_scores_gap = v_orig[1]
        
        # Predicted scores (placeholder/basic logic)
        pred_scores_list = predicted_scores if predicted_scores else []
        ppr_motif_scores = [] 
        if pred_scores_list:
             limit = min(motif_len, len(pred_scores_list))
             weighted = weighted_scores_vec[:limit] * np.array(pred_scores_list[:limit], dtype=np.float32)
             ppr_motif_scores = weighted.tolist()
             
        row = {
            "Gap_Index": front_size,
            "Gap_Size": gap_size,
            "Total_Score": score, # Weighted if use_ps_score
            "PPR_Binding_Score": score, # Consistency: Use the calculated score (weighted/unweighted)
            "Motif_Labels_GAP": ','.join(m_labels),
            "Motif_Types_GAP": ','.join(m_types), # Pass Types
            "Segment": sub_seq,
            "PPR_Binding_Scores_List": weighted_scores_vec.tolist(),
            "Motif_Confidence_Scores": conf_scores_gap,
            "PPR_Motif_Scores": ppr_motif_scores,
            "Predicted PPR Motif Scores List": pred_scores_list
        }
        results.append(row)
        
    return pd.DataFrame(results)

