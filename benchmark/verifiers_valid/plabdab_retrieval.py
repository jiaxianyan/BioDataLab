import os
from typing import Dict, Tuple

def plabdab_retrieval(d) -> bool:
    is_match = compare_fasta_files(
        output_fasta_path = f"{d}/antibody_seq_retrieval.fasta", 
    )
    return is_match

def compare_fasta_files(
    output_fasta_path: str, 
) -> bool:
    """
    Check if the output FASTA file exists.
    
    Args:
        output_fasta_path: Path to output FASTA file
    
    Returns:
        True if file exists, False otherwise
    """
    return os.path.exists(output_fasta_path)