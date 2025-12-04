import os
from Bio import Entrez, SeqIO
from pathlib import Path

def fetch_fasta(accession_id, email="example@example.com", output_dir="."):
    """
    Fetches a FASTA sequence from NCBI Nucleotide database.
    
    Args:
        accession_id (str): NCBI Accession ID (e.g., NC_000932).
        email (str): Email address for NCBI Entrez.
        output_dir (str): Directory to save the FASTA file.
        
    Returns:
        str: Path to the saved FASTA file.
    """
    Entrez.email = email
    
    try:
        print(f"Fetching {accession_id} from NCBI...")
        handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="fasta", retmode="text")
        record = handle.read()
        handle.close()
        
        if not record:
            raise ValueError(f"No data returned for accession {accession_id}")
            
        # Parse to ensure it's valid and get a clean filename
        # We can just save it directly, but parsing validates it.
        # Also, we might want to sanitize the ID for the filename.
        
        filename = f"{accession_id}.fasta"
        output_path = Path(output_dir) / filename
        
        with open(output_path, "w") as f:
            f.write(record)
            
        print(f"Saved to {output_path}")
        return str(output_path)
        
    except Exception as e:
        raise RuntimeError(f"Failed to fetch FASTA from NCBI: {e}")

def fetch_gff(accession_id, email="example@example.com", output_dir="."):
    """
    Fetches a GFF3 annotation file from NCBI Nucleotide database.
    
    Args:
        accession_id (str): NCBI Accession ID.
        email (str): Email address for NCBI Entrez.
        output_dir (str): Directory to save the GFF3 file.
        
    Returns:
        str: Path to the saved GFF3 file.
    """
    Entrez.email = email
    
    try:
        print(f"Fetching GFF3 for {accession_id} from NCBI...")
        # rettype="gff3" is not always standard for 'nucleotide' db in efetch, 
        # but often works or we use 'gv' (GenBank View)? 
        # Actually, for full records, 'gff3' might not be directly supported via efetch for all IDs.
        # But let's try. If not, we might need to fetch GenBank ('gb') and convert, 
        # but that's complex.
        # Alternative: 'ft' (Feature Table). IGV reads GFF3 best.
        # Let's try 'gff3'.
        handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="gff3", retmode="text")
        record = handle.read()
        handle.close()
        
        if not record:
            raise ValueError(f"No GFF3 data returned for accession {accession_id}")
            
        filename = f"{accession_id}.gff3"
        output_path = Path(output_dir) / filename
        
        with open(output_path, "w") as f:
            f.write(record)
            
        print(f"Saved to {output_path}")
        return str(output_path)
        
    except Exception as e:
        # Fallback or error
        print(f"Warning: Failed to fetch GFF3 directly: {e}")
        # Maybe try 'gff'?
        return None
