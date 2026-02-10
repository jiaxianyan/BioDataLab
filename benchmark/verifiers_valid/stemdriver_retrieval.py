import os

def stemdriver_retrieval(d):
    try:
        PRED_PATH = f"{d}/stemdriver_retrieval.json"
        return os.path.exists(PRED_PATH)
    except:
        return False