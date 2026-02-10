import os

def ddinter_annotate_1(d):
    try:
        # check if generated file exists
        file_path = f'{d}/ddinter_annotate_1.json'
        return os.path.exists(file_path)

    except Exception as e:
        print(f"Error checking file: {e}")
        return False