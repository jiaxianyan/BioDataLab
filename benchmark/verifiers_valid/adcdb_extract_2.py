import os

def adcdb_extract_2(d):
    try:
        PRED_PATH = f"{d}/adcdb_extract_2.json"
        return os.path.exists(PRED_PATH)
    except:
        return False