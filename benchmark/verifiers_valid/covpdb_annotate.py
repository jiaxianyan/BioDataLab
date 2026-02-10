import os

def covpdb_annotate(d):
    try:
        pred_path = os.path.join(d, "covpdb_annotate.txt")
        
        return os.path.exists(pred_path)
    
    except Exception as e:
        print(f"Error processing covpdb_annotate: {e}")
        return False