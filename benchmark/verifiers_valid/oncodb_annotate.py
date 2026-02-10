import os

def oncodb_annotate(d: str) -> bool:
    """
    Evaluate whether generated oncodb_annotate.csv exists.

    - Input:
        d: the run/output directory root (where <r></r>/oncodb_annotate.csv is saved)

    - Logic:
        * Check if the generated file exists using os.path.exists()
        * Return True if exists, False otherwise
    """
    try:
        pred_path = os.path.join(d, "oncodb_annotate.csv")
        
        if os.path.exists(pred_path):
            return True
        else:
            print(f"Missing prediction file: {pred_path}")
            return False
    
    except Exception as e:
        print(f"Error evaluating oncodb_annotate: {e}")
        return False