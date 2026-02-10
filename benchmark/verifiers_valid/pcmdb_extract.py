import os

def pcmdb_extract(d):
    try:
        PRED_PATH = f"{d}/pcmdb_extract.csv"
        return os.path.exists(PRED_PATH)
    except:
        return False