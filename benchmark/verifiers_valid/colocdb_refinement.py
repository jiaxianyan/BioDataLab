import os

def colocdb_refinement(d: str) -> bool:
    """
    Check whether `{d}/colocdb_refinement.tsv` exists.
    
    Returns True if the file exists, False otherwise.
    """
    pred_path = f"{d}/colocdb_refinement.tsv"
    return os.path.exists(pred_path)