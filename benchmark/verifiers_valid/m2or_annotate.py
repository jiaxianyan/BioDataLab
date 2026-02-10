import os

def m2or_annotate(d):
    try:
        # check if predicted result file exists
        return os.path.exists(f'{d}/m2or_annotate.json')
    
    except Exception as e:
        print(f"Error checking file: {e}")
        return False