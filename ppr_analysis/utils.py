"""
Utility functions for PPR analysis.
"""
import pandas as pd

def infer_additional_info(file_name):
    """
    Infers analysis parameters from the filename.
    
    Args:
        file_name (str): The name of the file.
        
    Returns:
        tuple: (gap_info, motif_info, binding_source)
    """
    gap_info = "GAP not considered" if "gap0" in file_name else "GAP considered"
    motif_info = "Motif score not included" if "no_ps" in file_name else "Motif score included"
    binding_source = "Yan" if "yan" in file_name.lower() else "Kobayashi"
    return gap_info, motif_info, binding_source

def get_top_n_per_ppr(df, column_name, top_n=10):
    """
    Returns the top N rows per PPR_Name based on a specific column.
    
    Args:
        df (pd.DataFrame): Input DataFrame.
        column_name (str): Column to sort by.
        top_n (int): Number of top rows to return.
        
    Returns:
        pd.DataFrame: Filtered DataFrame.
    """
    top_n_df = df.groupby('PPR_Name').apply(lambda x: x.nlargest(top_n, column_name)).reset_index(drop=True)
    return top_n_df
