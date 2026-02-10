import os

def themarker_extract_2(d):
    try:
        PRED_PATH = f"{d}/themarker_extract_2.csv"
        return os.path.exists(PRED_PATH)
    except:
        return False