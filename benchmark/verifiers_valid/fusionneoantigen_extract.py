import os

def fusionneoantigen_extract(d):
    try:
        PRED_PATH = f"{d}/fusionneoantigen_extract.json"
        return os.path.exists(PRED_PATH)
    except:
        return False