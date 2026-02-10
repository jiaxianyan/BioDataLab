import os

def bioka_extract(d):
    try:
        pred_path = f"{d}/bioka_extract.csv"
        return os.path.exists(pred_path)
    except Exception as e:
        print(f"Error checking file: {e}")
        return False