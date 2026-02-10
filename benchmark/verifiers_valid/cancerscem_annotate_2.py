import os

def cancerscem_annotate_2(d):
    try:
        pred_file = f"{d}/cancerscem_annotate_2.txt"
        return os.path.exists(pred_file)
    except Exception as e:
        print(f"Error processing cancerscem_annotate_2 results: {e}")
        return False