import os

def scqtlbase_retrieval(d):
    try:
        PRED_PATH = f"{d}/scqtlbase_retrieval.json"
        return os.path.exists(PRED_PATH)
    except:
        return False