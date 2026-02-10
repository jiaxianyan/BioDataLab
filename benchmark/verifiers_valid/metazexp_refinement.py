import os

def metazexp_refinement(d):
    try:
        # Check if generated file exists
        return os.path.exists(f'{d}/metazexp_refinement.json')
    
    except Exception as e:
        print(f"Error processing metazexp_refinement evaluation: {e}")
        return False