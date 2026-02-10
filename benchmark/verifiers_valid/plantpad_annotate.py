import os

def plantpad_annotate(d):
    pred_path = os.path.join(d, "plantpad_annotate.csv")
    return os.path.exists(pred_path)