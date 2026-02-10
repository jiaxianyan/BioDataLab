import os

def a3d_modb_retrieval(d):
    try:
        # check if predicted result file exists
        return os.path.exists(f'{d}/a3d_modb_retrieval.json')
    except Exception as e:
        print(f"Error checking file: {e}")
        return False