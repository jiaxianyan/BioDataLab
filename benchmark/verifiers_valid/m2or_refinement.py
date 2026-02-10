import os

def m2or_refinement(d):
    """
    Check whether the generated file exists.
    
    Args:
        d (str): directory path that contains the model output file
    
    Returns:
        bool: True if file exists, False otherwise
    """
    # path to predicted result
    pred_path = f"{d}/m2or_refinement.txt"
    
    return os.path.exists(pred_path)