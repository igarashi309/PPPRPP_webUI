"""
Input/Output operations for PPR analysis.
"""
import pickle
import csv
import os

def read_scores_from_file(filename):
    """
    Reads scores from a text file.
    
    Args:
        filename (str): Path to the file.
        
    Returns:
        list: List of float scores.
    """
    if not os.path.exists(filename):
        # Fallback or error handling
        print(f"Warning: Score file {filename} not found. Returning empty list.")
        return []
        
    with open(filename, 'r') as f:
        scores = f.readlines()
    scores = [float(score.strip()) for score in scores]
    return scores

def load_ppr_data_from_pickle(pickle_file):
    """
    Loads PPR data from a pickle file.
    Mimics the notebook logic: extracts PPR_CODE and PS_SCORE from annotations.
    
    Args:
        pickle_file (str or file-like): Path to the pickle file or file object.
        
    Returns:
        dict: {'motif_labels': [], 'ppr_scores': []}
    """
    if isinstance(pickle_file, str):
        with open(pickle_file, 'rb') as f:
            data = pickle.load(f)
    else:
        # Assume it's a file-like object (e.g. from Streamlit)
        data = pickle.load(pickle_file)
        
    motif_labels = []
    ppr_scores = []
    
    # Assuming data is a list of objects with .annotations attribute
    for p in data:
        if "PPR_CODE" in p.annotations and "PS_SCORE" in p.annotations:
            code = p.annotations["PPR_CODE"]
            score = p.annotations["PS_SCORE"]
            
            # Split code into 2-char motifs
            motifs = [code[i:i+2] for i in range(0, len(code), 2)]
            motif_labels.extend(motifs)
            
            # Replicate score for each motif
            ppr_scores.extend([score] * len(motifs))
            
    return {'motif_labels': motif_labels, 'ppr_scores': ppr_scores}

def save_scores_to_csv(all_scores, output_dir="."):
    """
    Saves aggregated scores to CSV files using standard csv module.
    
    Args:
        all_scores (dict): Dictionary of scores keyed by analysis option.
        output_dir (str): Directory to save CSVs.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    headers = [
        "Position", "PPR_Binding_Score", "PPR_Name", "Motif_Labels", 
        "Motif_Labels_GAP", "Gap_Size", "Segment", "PPR_Motif_Scores", 
        "Total_Score", "PPR_Binding_Scores_List", "Predicted PPR Motif Scores List"
    ]

    for key, data in all_scores.items():
        output_path = os.path.join(output_dir, f"{key}.csv")
        
        # data is now a dict of columns
        # We need to iterate over rows
        if not data or not data["Position_Relative"]:
            continue
            
        num_rows = len(data["Position_Relative"])
        
        with open(output_path, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(headers)
            
            # Columns in order of headers
            # Headers: Position, PPR_Binding_Score, PPR_Name, Motif_Labels, ...
            # Map headers to data keys
            # Note: Headers in CSV are slightly different from keys in dict?
            # Dict keys: Position_Relative, PPR_Binding_Score, PPR_Name, Motif_Labels, Motif_Labels_GAP, Gap_Size, Segment, PPR_Motif_Scores, Total_Score, PPR_Binding_Scores_List, Predicted PPR Motif Scores List
            # CSV Headers: Position, PPR_Binding_Score, PPR_Name, Motif_Labels, Motif_Labels_GAP, Gap_Size, Segment, PPR_Motif_Scores, Total_Score, PPR_Binding_Scores_List, Predicted PPR Motif Scores List
            
            # Mapping
            col_map = [
                "Position_Relative", "PPR_Binding_Score", "PPR_Name", "Motif_Labels", 
                "Motif_Labels_GAP", "Gap_Size", "Segment", "PPR_Motif_Scores", 
                "Total_Score", "PPR_Binding_Scores_List", "Predicted PPR Motif Scores List"
            ]
            
            # Prepare columns for iteration
            cols = [data[k] for k in col_map]
            
            # Iterate
            for i in range(num_rows):
                row = [c[i] for c in cols]
                
                # Format lists as strings
                # Indices: 7 (PPR_Motif_Scores), 9 (PPR_Binding_Scores_List), 10 (Predicted...)
                row[7] = ','.join(map(str, row[7]))
                row[9] = ','.join(map(str, row[9]))
                row[10] = ','.join(map(str, row[10]))
                
                writer.writerow(row)
                
        print(f"Saved {output_path}")
