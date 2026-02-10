import os

def circmine_retrieval(d):
    try:
        PRED_PATH = f"{d}/circmine_retrieval.json"
        return os.path.exists(PRED_PATH)
    except:
        return False