import os

def covpdb_integration(d,
             pred_name="covpdb_integration_complex.csv",
             gold_path="benchmark/gold_results/covpdb_integration_complex.csv",
             rtol=1e-3,
             atol=1e-2,
             affinity_rtol=1e-2,
             affinity_atol=1e-2):
    """
    Check whether predicted csv file exists.

    Args:
        d (str): directory containing prediction csv (d/covpdb_integration_complex.csv)
        pred_name (str): prediction filename
        gold_path (str): ground truth csv path (unused)
        rtol/atol: tolerance parameters (unused)
        affinity_rtol/affinity_atol: tolerance parameters (unused)

    Returns:
        bool: True if prediction file exists, False otherwise
    """
    try:
        pred_path = os.path.join(d, pred_name)
        if os.path.exists(pred_path):
            return True
        else:
            print(f"[EVAL] Pred file not found: {pred_path}")
            return False
    except Exception as e:
        print(f"[EVAL] Error checking file: {e}")
        return False