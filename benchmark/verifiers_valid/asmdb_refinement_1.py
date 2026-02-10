import os

def asmdb_refinement_1(d):
    try:
        pred_path = f"{d}/asmdb_refinement_1.bed"
        
        return os.path.exists(pred_path)

    except Exception as e:
        print(f"Error evaluating asmdb_refinement_1: {e}")
        return False