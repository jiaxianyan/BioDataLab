import os

def crost_retrieval(d):
    try:
        PRED_PATH = f"{d}/crost_retrieval.json"
        return os.path.exists(PRED_PATH)
    except:
        return False