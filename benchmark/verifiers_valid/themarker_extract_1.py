import os

def themarker_extract_1(d):
    try:
        PRED_PATH = f"{d}/themarker_extract_1.json"
        return os.path.exists(PRED_PATH)
    except:
        return False