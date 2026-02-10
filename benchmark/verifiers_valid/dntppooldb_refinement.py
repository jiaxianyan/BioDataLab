import os


def dntppooldb_refinement(d: str) -> bool:
    """
    Verifier for dntppooldb_refinement.

    Checks if the generated file exists.
    """
    try:
        pred_path = os.path.join(d, "dntppooldb_refinement.csv")
        
        if os.path.exists(pred_path):
            return True
        else:
            print(f"Missing prediction file: {pred_path}")
            return False

    except Exception as e:
        print(f"Error: {e}")
        return False