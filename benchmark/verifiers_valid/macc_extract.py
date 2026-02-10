import os

def macc_extract(d):
    try:
        PRED_PATH = f"{d}/macc_extract.csv"
        return os.path.exists(PRED_PATH)
    except:
        return False