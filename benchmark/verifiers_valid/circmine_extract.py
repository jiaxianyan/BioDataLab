import os

def circmine_extract(d):
    try:
        PRED_PATH = f"{d}/circmine_extract.csv"
        return os.path.exists(PRED_PATH)
    except:
        return False