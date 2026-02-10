import os

def scovid_refinement(d: str) -> bool:
    """
    Evaluate whether the predicted qc_passed_barcodes.csv exists.

    Pred file:  {d}/qc_passed_barcodes.csv

    Check:
      - Only verify if the file exists
    """
    try:
        pred_path = os.path.join(d, "qc_passed_barcodes.csv")
        
        if os.path.exists(pred_path):
            return True
        else:
            print(f"[scov2_md_annotate] Missing prediction file: {pred_path}")
            return False

    except Exception as e:
        print(f"[scov2_md_annotate] Error: {e}")
        return False