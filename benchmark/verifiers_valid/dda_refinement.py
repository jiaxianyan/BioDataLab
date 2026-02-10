import os

def dda_refinement(d):
    try:
        pred_path = f"{d}/dda_refinement.fastq.gz"
        return os.path.exists(pred_path)
    except Exception as e:
        print(f"Error checking file: {e}")
        return False