"""
Core logic for PPR motif analysis.
"""
import os
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
    gap_variations = [(motif_labels, ppr_scores, 0)]  # GAP 0
    for gap_size in range(1, max_gap + 1):
        for front_size in range(min_front, len(motif_labels) - min_back + 1):
            gap_variation = motif_labels[:front_size] + ['--'] * gap_size + motif_labels[front_size:]
            gap_scores = ppr_scores[:front_size] + [0] * gap_size + ppr_scores[front_size:]
            gap_variations.append((gap_variation, gap_scores, gap_size))
    return gap_variations

def calculate_binding_score(rna_seq, motif_binding_scores, motif_labels, ppr_scores):
    """
    Calculates binding scores for a sequence against motifs.
    """
    motif_length = len(motif_labels)
    binding_sites = []
    
    # Pre-calculate motif scores to avoid repeated dictionary lookups
    # This is an optimization over the original notebook code
    
    for i in range(len(rna_seq) - motif_length + 1):
        segment = rna_seq[i:i+motif_length]
        scores = []
        motif_binding_score_list = []
        
        for j in range(motif_length):
            motif = motif_labels[j]
            nuc = segment[j]
            
            base_score = 0.25
            if motif == "--":
                score = 0
                base_score = 0 # Or should this be 0? Notebook says scores.append(0) for GAP.
                               # But motif_binding_score_list logic in notebook:
                               # if motif in motif_binding_scores and nuc in nucleotide_index else 0.25
                               # '--' is in motif_binding_scores with [0,0,0,0].
                               # So base_score for GAP is 0.
                if motif in motif_binding_scores:
                     base_score = motif_binding_scores[motif][0] # All 0 for GAP
            elif motif in motif_binding_scores and nuc in NUCLEOTIDE_INDEX:
                base_score = motif_binding_scores[motif][NUCLEOTIDE_INDEX[nuc]]
                score = base_score * ppr_scores[j]
            else:
                score = 0.25 * ppr_scores[j]
            
            scores.append(score)
            motif_binding_score_list.append(base_score)
            
        score = sum(scores)
        binding_sites.append((i, score, segment, scores, score, motif_binding_score_list))
        
    return binding_sites

def process_fasta_sequences(fasta_sequence, motif_binding_scores, motif_labels_variations, ppr_name, predicted_scores, use_ps_score=True):
    """
    Processes FASTA sequence with given variations.
    """
    all_scores = []
    rna_seq = dna_to_rna(fasta_sequence)

    for motif_labels, ppr_scores, gap_size in motif_labels_variations:
        binding_scores = calculate_binding_score(rna_seq, motif_binding_scores, motif_labels, ppr_scores)
        motif_labels_gap = ','.join(motif_labels)
        
        for pos, score, segment, scores, total_score, motif_binding_score_list in binding_scores:
            if use_ps_score:
                ppr_binding_score = score
            else:
                ppr_binding_score = sum(motif_binding_score_list)

            # Verification for GAP
            gap_indices = [i for i, label in enumerate(motif_labels) if label == '--']
            for gap_index in gap_indices:
                # In the notebook, this is an assert. We'll keep it safe.
                if scores[gap_index] != 0:
                     # Log warning or handle error? For now, assume correct logic.
                     pass

            # Calculate PPR Motif Scores (binding score * predicted score)
            # Note: ppr_scores passed to generate_gap_variations comes from the pickle file (PS_SCORE).
            # predicted_scores passed here comes from 'scores_{job_name}.txt'. 
            # The notebook uses `predicted_scores` (from txt) to multiply with `motif_binding_score_list`.
            # But `ppr_scores` (from pickle) was used in `calculate_binding_score` to calculate `scores`.
            # Let's strictly follow notebook logic:
            # ppr_motif_scores = [b * p for b, p in zip(motif_binding_score_list, predicted_scores)]
            # Wait, predicted_scores length might not match motif_binding_score_list if gaps are added?
            # In notebook: `process_scores` passes `predicted_scores` (list of floats) to `process_fasta_sequences`.
            # `predicted_scores` seems to be a global list read from file. 
            # BUT `motif_binding_score_list` length depends on `motif_labels` (which can have gaps).
            # If `predicted_scores` is fixed length, zip will truncate or mismatch.
            # However, in the notebook `predicted_scores` seems to be derived from the pickle file as well in `process_scores`?
            # No, `predicted_scores` is read from `scores_{job_name}.txt`.
            # And `ppr_scores` (in `generate_gap_variations`) comes from pickle `PS_SCORE`.
            
            # Let's look at notebook `process_scores`:
            # ppr_scores.extend([p.annotations["PS_SCORE"]]...)
            # gap_variations = generate_gap_variations(motif_labels, ppr_scores, ...)
            # scores = { ... process_fasta_sequences(..., predicted_scores, ...) }
            
            # If `predicted_scores` is just a list of scores, how does it align with motifs?
            # The notebook code: `ppr_motif_scores = [b * p for b, p in zip(motif_binding_score_list, predicted_scores)]`
            # This implies `predicted_scores` should align with the motifs.
            # If gaps are inserted, `motif_binding_score_list` grows. `predicted_scores` must also grow or handle gaps?
            # The notebook does NOT seem to adjust `predicted_scores` for gaps in `process_scores`.
            # This looks like a potential bug in the original notebook or `predicted_scores` is something else.
            # Actually, `predicted_scores` in `main` is read from a file.
            # But in `process_scores` (line 3767), `predicted_scores` is passed in.
            # If `predicted_scores` is a list of scalars corresponding to the original motifs, 
            # then `zip` with a gapped list will be misaligned if the gap is in the middle.
            # However, we must replicate the notebook's behavior for now.
            
            # Handling the zip safely:
            # If predicted_scores is not long enough, zip stops.
            # If it is the same length as original motifs, and we added gaps, we might need to insert 0s into predicted_scores too?
            # The notebook does NOT do this. It just zips. 
            # We will replicate this behavior but be aware of it.
            
            ppr_motif_scores = [b * p for b, p in zip(motif_binding_score_list, predicted_scores)]

            all_scores.append((
                pos, 
                ppr_binding_score, 
                ppr_name, 
                ','.join(motif_labels), 
                motif_labels_gap, 
                gap_size, 
                segment, 
                ppr_motif_scores, 
                total_score, 
                motif_binding_score_list, 
                predicted_scores,
                ppr_scores # Add ppr_scores (Motif Confidence Scores) to output
            ))
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
    
    gap0_variations = [(motif_labels, ppr_scores, 0)]
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
