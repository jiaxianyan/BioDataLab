import os

def adcdb_extract_1(d):
    try:
        PRED_PATH = f"{d}/adcdb_extract_1.csv"
        return os.path.exists(PRED_PATH)
    except:
        return False