import os

def ravar_extract(d):
    try:
        PRED_PATH = f"{d}/ravar_extract.json"
        return os.path.exists(PRED_PATH)
    except:
        return False