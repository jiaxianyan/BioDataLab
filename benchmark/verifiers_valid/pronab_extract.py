import os

def pronab_extract(d):
    try:
        PRED_PATH = f"{d}/pronab_extract.csv"
        return os.path.exists(PRED_PATH)
    except:
        return False