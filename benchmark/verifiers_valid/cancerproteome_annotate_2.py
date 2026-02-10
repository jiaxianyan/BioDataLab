import os

def cancerproteome_annotate_2(d):
    """
    Evaluate whether the generated cancerproteome_annotate_2.txt exists.
    - Input: d (str) directory path that contains the predicted result file.
    - Prediction file: <d>/cancerproteome_annotate_2.txt
    - Returns: True if file exists, False otherwise.
    """
    pred_path = os.path.join(d, "cancerproteome_annotate_2.txt")
    
    return os.path.exists(pred_path)