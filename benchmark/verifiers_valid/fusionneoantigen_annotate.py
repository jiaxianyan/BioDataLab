import os

def fusionneoantigen_annotate(d):
    """
    Check whether the predicted fusion neoantigen peptide file exists.
    
    Logic:
    - Check if the generated file exists using os.path.exists()
    - Return True if exists, False otherwise
    """
    try:
        pred_fasta = f"{d}/fusionneoantigen_annotate.fasta"
        return os.path.exists(pred_fasta)
    except Exception as e:
        print(f"Error evaluating fusionneoantigen_annotate: {e}")
        return False