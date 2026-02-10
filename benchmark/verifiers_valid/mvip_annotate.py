import os

def mvip_annotate(d):
    try:
        # check if generated file exists
        return os.path.exists(f'{d}/mvip_annotate.tsv')
    
    except Exception as e:
        print(f"Error processing mvip_annotate results: {e}")
        return False