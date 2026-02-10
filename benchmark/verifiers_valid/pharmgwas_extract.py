import os

def pharmgwas_extract(d):
    try:
        PRED_PATH = f"{d}/pharmgwas_extract.csv"
        return os.path.exists(PRED_PATH)
    except:
        return False