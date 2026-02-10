import os

def themarker_annotate(d):
    try:
        # Check if predicted results file exists
        pred_file = f'{d}/themarker_annotate.json'
        return os.path.exists(pred_file)

    except Exception as e:
        print(f"Error checking file: {e}")
        return False