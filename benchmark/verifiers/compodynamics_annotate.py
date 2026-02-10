import os
import math

def compodynamics_annotate(d, rel_tol=1e-6, abs_tol=1e-9):
    """
    Evaluate whether the predicted genome-wide weighted GC content equals the ground truth.

    Expected prediction file:
      {d}/compodynamics_annotate.txt

    Expected gold file:
      benchmark/gold_results/compodynamics_annotate.txt

    Comparison:
      float values compared with math.isclose (relative + absolute tolerance).
    """
    def _read_float(path: str) -> float:
        with open(path, "r", encoding="utf-8") as f:
            s = f.read().strip()
        # be tolerant to minor formatting (spaces, trailing newline, accidental commas)
        s = s.replace(",", "").strip()
        if not s:
            raise ValueError(f"Empty result file: {path}")
        # take the first token that looks like a number
        token = s.split()[0]
        return float(token)

    try:
        pred_path = os.path.join(d, "compodynamics_annotate.txt")
        gold_path = os.path.join("benchmark", "gold_results", "compodynamics_annotate.txt")

        pred_val = _read_float(pred_path)
        gold_val = _read_float(gold_path)

        return math.isclose(pred_val, gold_val, rel_tol=rel_tol, abs_tol=abs_tol)

    except FileNotFoundError as e:
        print(f"Missing file: {e}")
        return False
    except Exception as e:
        print(f"Error evaluating compodynamics_annotate: {e}")
        return False
