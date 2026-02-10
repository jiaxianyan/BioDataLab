import os

def tf_marker_retrieval(d):
    try:
        PRED_PATH = f"{d}/tf_marker_retrieval.json"
        return os.path.exists(PRED_PATH)
    except:
        return False