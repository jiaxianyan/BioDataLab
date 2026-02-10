import os

def ddinter_integration_1(d):
    try:
        pred_path = f"{d}/ddinter_integration_1.csv"
        
        return os.path.exists(pred_path)
    
    except Exception as e:
        print(f"Error checking file: {e}")
        return False