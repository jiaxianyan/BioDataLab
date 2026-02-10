import os

def amdb_retrieval(d):
    return os.path.exists(f'{d}/amdb_retrieval.json')