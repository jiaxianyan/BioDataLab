import os

def compodynamics_integration(d):
    pred_fna_path = f"{d}/compodynamics_integration.fna"
    return os.path.exists(pred_fna_path)