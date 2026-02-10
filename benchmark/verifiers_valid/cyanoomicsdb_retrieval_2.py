import os

def cyanoomicsdb_retrieval_2(d):
    PRED_PATH = f"{d}/cyanoomicsdb_retrieval_2.csv"
    return os.path.exists(PRED_PATH)