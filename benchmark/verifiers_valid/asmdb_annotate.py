import os

def asmdb_annotate(d):
    try:
        pred_path = f"{d}/asmdb_annotate.csv"
        return os.path.exists(pred_path)
    except Exception as e:
        print(f"Error checking file: {e}")
        return False