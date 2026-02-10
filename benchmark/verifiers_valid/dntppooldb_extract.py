import os

def dntppooldb_extract(d):
    PRED_PATH = f"{d}/dntppooldb_extract.csv"
    return os.path.exists(PRED_PATH)