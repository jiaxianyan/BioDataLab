import os

def cyanoomicsdb_annotate_1(d: str) -> bool:
    """
    Evaluate whether the generated file exists.

    Pred path:  {d}/cyanoomicsdb_annotate_1.tsv
    """
    pred_path = os.path.join(d, "cyanoomicsdb_annotate_1.tsv")
    return os.path.exists(pred_path)