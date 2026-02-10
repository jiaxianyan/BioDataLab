import os

def metazexp_annotate(d):
    """
    Checks if the generated file exists.
    
    Args:
        d (str): Directory containing the predicted 'metazexp_annotate.tsv'.
        
    Returns:
        bool: True if the file exists, False otherwise.
    """
    PRED_FILENAME = 'metazexp_annotate.tsv'
    pred_path = os.path.join(d, PRED_FILENAME)
    
    return os.path.exists(pred_path)