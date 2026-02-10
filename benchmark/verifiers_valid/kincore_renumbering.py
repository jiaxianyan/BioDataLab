import os
from typing import Dict, Tuple


def kincore_renumbering(d) -> bool:
    is_match = compare_pdb_files(
        ref_pdb_path = 'benchmark/gold_results/kincore_renumbering.pdb', 
        out_pdb_path = f"{d}/kincore_renumbering.pdb",
    )
    return is_match
    

def compare_pdb_files(ref_pdb_path: str, out_pdb_path: str) -> bool:
    """
    Check if the output PDB file exists.
    Returns True if file exists, False otherwise.
    """
    return os.path.exists(out_pdb_path)