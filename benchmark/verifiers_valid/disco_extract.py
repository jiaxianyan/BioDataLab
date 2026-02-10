import os

def disco_extract(d):
    try:
        PRED_PATH = f"{d}/disco_extract.csv"
        return os.path.exists(PRED_PATH)
    except:
        return False