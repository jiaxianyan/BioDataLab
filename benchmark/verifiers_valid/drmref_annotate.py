import os

def drmref_annotate(d):
    """
    Check whether `<d>/drmref_annotate.csv` exists.
    
    Returns True if the file exists, False otherwise.
    """
    try:
        pred_path = f"{d}/drmref_annotate.csv"
        return os.path.exists(pred_path)
    except Exception as e:
        print(f"Error checking file: {e}")
        return False