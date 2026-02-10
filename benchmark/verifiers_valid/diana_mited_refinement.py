import os

def diana_mited_refinement(d, threshold=0.5):
    """
    Check whether the predicted result file exists.

    Args:
        d (str): directory containing the predicted result file
        threshold (float): allowed absolute difference in percentage points (unused)

    Returns:
        bool: True if file exists, False otherwise
    """
    pred_path = f"{d}/diana_mited_refinement.txt"
    return os.path.exists(pred_path)