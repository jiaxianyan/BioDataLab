import os

def scapaatlas_annotate(d: str) -> bool:
    """
    Evaluate whether the generated file exists.
    - Input:
        d: directory containing the prediction file `scapaatlas_annotate.bed`
    - Check:
        pred: {d}/scapaatlas_annotate.bed
    - Returns:
        True if file exists, False otherwise
    """
    PRED_NAME = "scapaatlas_annotate.bed"
    pred_path = os.path.join(d, PRED_NAME)
    
    return os.path.exists(pred_path)