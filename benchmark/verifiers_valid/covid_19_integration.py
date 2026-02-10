import os

def covid_19_integration(d):
    try:
        pred_path = f"{d}/covid_19_integration.csv"
        return os.path.exists(pred_path)
    except Exception as e:
        print(f"Error checking file: {e}")
        return False