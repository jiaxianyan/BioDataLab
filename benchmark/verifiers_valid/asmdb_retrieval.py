import os

def asmdb_retrieval(d):
    try:
        PRED_PATH = f"{d}/asmdb_retrieval.json"
        return os.path.exists(PRED_PATH)
    except:
        return False