"""
NCBI Integration module for fetching GenBank records.
"""
from Bio import Entrez, SeqIO
from io import StringIO

def fetch_genbank_record(accession_id, email="example@email.com"):
    """
    Fetches a GenBank record by accession ID.
    
    Args:
        accession_id (str): GenBank Accession ID (e.g., AP012527).
        email (str): Email for Entrez identification.
        
    Returns:
        Bio.SeqRecord.SeqRecord: The fetched record, or None if failed.
    """
    Entrez.email = email
    try:
        # Fetch GenBank file
        handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        return record
    except Exception as e:
        print(f"Error fetching GenBank record {accession_id}: {e}")
        return None

def extract_features(record):
    """
    Extracts relevant features (CDS, genes) from a SeqRecord.
    
    Args:
        record (Bio.SeqRecord.SeqRecord): The record.
        
    Returns:
        list: List of dicts {'type': str, 'start': int, 'end': int, 'strand': int, 'name': str}
    """
    features = []
    for feature in record.features:
        if feature.type in ['CDS', 'gene', 'rRNA', 'tRNA']:
            start = int(feature.location.start)
            end = int(feature.location.end)
            strand = feature.location.strand
            
            name = feature.qualifiers.get('gene', [''])[0]
            if not name:
                name = feature.qualifiers.get('product', [''])[0]
                
            features.append({
                'type': feature.type,
                'start': start,
                'end': end,
                'strand': strand, # 1 for forward, -1 for reverse
                'name': name
            })
    return features
