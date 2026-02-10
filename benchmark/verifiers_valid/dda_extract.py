import os

def dda_extract(d):
    try:
        PRED_PATH = f"{d}/dda_extract.csv"
        return os.path.exists(PRED_PATH)
    except:
        return False