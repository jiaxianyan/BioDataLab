import os

def pncshub_annotate(d):
    try:
        pred_path = os.path.join(d, "pncshub_annotate.txt")
        return os.path.exists(pred_path)
    except Exception as e:
        print(f"Error processing pncshub_annotate text files: {e}")
        return False