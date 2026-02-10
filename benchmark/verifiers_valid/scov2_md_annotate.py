import os

def scov2_md_annotate(d, abs_tol=1e-3, rel_tol=1e-2):
    """
    Evaluate whether predicted file exists.

    Expected prediction path:
      {d}/rmsf_ca.json

    Logic:
      Check if the prediction file exists using os.path.exists()
    """
    pred_path = os.path.join(d, "rmsf_ca.json")
    return os.path.exists(pred_path)