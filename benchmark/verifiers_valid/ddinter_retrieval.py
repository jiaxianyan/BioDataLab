import os

def ddinter_retrieval(d):
    try:
        PRED_PATH = f"{d}/ddinter_retrieval.json"
        return os.path.exists(PRED_PATH)
    except:
        return False