import os

def ddinter_annotate_2(d):
    try:
        # check if generated file exists
        return os.path.exists(f'{d}/ddinter_annotate_2.json')
    
    except Exception as e:
        print(f"Error checking file: {e}")
        return False