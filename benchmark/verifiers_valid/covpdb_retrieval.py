import os

def covpdb_retrieval(d):
    try:
        file_path = f'{d}/cov_pdb_retrieval.json'
        return os.path.exists(file_path)
        
    except Exception as e:
        print(f"Error checking file: {e}")
        return False