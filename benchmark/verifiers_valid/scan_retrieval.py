import os

def scan_retrieval(d):
    try:
        PRED_PATH = f"{d}/scan_retrieval.json"
        return os.path.exists(PRED_PATH)
    except:
        return False