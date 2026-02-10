import os

def gpedit_refinement(d: str) -> bool:
    """
    Verifier for GPEdit task:
    Check if the generated file exists.
    
    Returns:
    - True if <d>/gpedit_refinement.tsv exists
    - False otherwise
    """
    pred_path = f"{d}/gpedit_refinement.tsv"
    return os.path.exists(pred_path)