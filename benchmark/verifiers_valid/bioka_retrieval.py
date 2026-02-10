import os

def bioka_retrieval(d):
    try:
        PRED_PATH = f"{d}/bioka_retrieval.json"
        return os.path.exists(PRED_PATH)
    except:
        return False