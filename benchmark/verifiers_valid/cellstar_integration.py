import os

def cellstar_integration(d):
    try:
        pred_path = f"{d}/cellstar_integration.csv"
        return os.path.exists(pred_path)
    except Exception as e:
        print(f"Error evaluating cellstar_integration: {e}")
        return False