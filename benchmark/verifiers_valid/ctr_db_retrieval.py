import os

def ctr_db_retrieval(d):
    try:
        PRED_PATH = f"{d}/ctr_db_retrieval.json"
        return os.path.exists(PRED_PATH)
    except:
        return False