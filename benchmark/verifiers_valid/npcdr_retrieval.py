import os

def npcdr_retrieval(d):
    try:
        PRED_PATH = f"{d}/npcdr_retrieval.json"
        return os.path.exists(PRED_PATH)
    except:
        return False