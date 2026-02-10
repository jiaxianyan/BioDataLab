import json
import math
import os

def scov2_md_annotate(d, abs_tol=1e-3, rel_tol=1e-2):
    """
    Evaluate whether predicted RMSF-CA JSON equals ground truth.

    Expected prediction path:
      {d}/rmsf_ca.json

    Ground truth path (tries multiple common locations):
      benchmark/gold_results/scov2_md_annotate/rmsf_ca.json
      benchmark/gold_results/scov2_md_annotate.json
      benchmark/gold_results/rmsf_ca.json

    Logic:
      1) Load pred + gold
      2) Normalize keys to stringified 0-based ints
      3) Check total count and key set
      4) Compare each value with tolerance (abs + rel)
    """
    try:
        pred_path = os.path.join(d, "rmsf_ca.json")
        if not os.path.exists(pred_path):
            print(f"[scov2_md_annotate] Missing prediction file: {pred_path}")
            return False

        # Try a few likely gold locations to be robust
        gold_candidates = [
            os.path.join("benchmark", "gold_results", "scov2_md_annotate", "rmsf_ca.json"),
            os.path.join("benchmark", "gold_results", "scov2_md_annotate.json"),
            os.path.join("benchmark", "gold_results", "rmsf_ca.json"),
        ]
        gold_path = next((p for p in gold_candidates if os.path.exists(p)), None)
        if gold_path is None:
            print("[scov2_md_annotate] Missing ground truth file. Tried:\n  - " + "\n  - ".join(gold_candidates))
            return False

        with open(pred_path, "r", encoding="utf-8") as f:
            pred = json.load(f)
        with open(gold_path, "r", encoding="utf-8") as f:
            gold = json.load(f)

        if not isinstance(pred, dict) or not isinstance(gold, dict):
            print("[scov2_md_annotate] Prediction and gold must both be JSON objects (dict).")
            return False

        def normalize(obj):
            norm = {}
            for k, v in obj.items():
                # normalize key -> "int"
                try:
                    ki = int(k) if not isinstance(k, int) else k
                except Exception:
                    print(f"[scov2_md_annotate] Non-integer residue index key found: {k!r}")
                    return None
                ks = str(ki)

                # normalize value -> float
                try:
                    vf = float(v)
                except Exception:
                    print(f"[scov2_md_annotate] Non-numeric RMSF value for residue {ks}: {v!r}")
                    return None

                # disallow NaN/inf
                if not math.isfinite(vf):
                    print(f"[scov2_md_annotate] Non-finite RMSF value for residue {ks}: {vf}")
                    return None

                norm[ks] = vf
            return norm

        pred_n = normalize(pred)
        gold_n = normalize(gold)
        if pred_n is None or gold_n is None:
            return False

        # 1) quantity check
        if len(pred_n) != len(gold_n):
            print(f"[scov2_md_annotate] Count mismatch: pred={len(pred_n)} gold={len(gold_n)}")
            # still check key differences to help debugging
            pred_keys = set(pred_n.keys())
            gold_keys = set(gold_n.keys())
            missing = sorted(gold_keys - pred_keys)[:10]
            extra = sorted(pred_keys - gold_keys)[:10]
            if missing:
                print(f"[scov2_md_annotate] Missing keys (showing up to 10): {missing}")
            if extra:
                print(f"[scov2_md_annotate] Extra keys (showing up to 10): {extra}")
            return False

        # 2) key set check
        pred_keys = set(pred_n.keys())
        gold_keys = set(gold_n.keys())
        if pred_keys != gold_keys:
            missing = sorted(gold_keys - pred_keys)[:10]
            extra = sorted(pred_keys - gold_keys)[:10]
            print("[scov2_md_annotate] Key set mismatch.")
            if missing:
                print(f"[scov2_md_annotate] Missing keys (showing up to 10): {missing}")
            if extra:
                print(f"[scov2_md_annotate] Extra keys (showing up to 10): {extra}")
            return False

        # 3) numeric comparison (approx allowed)
        # Use a stable iteration order for reproducible error messages
        for k in sorted(gold_keys, key=lambda x: int(x)):
            pv = pred_n[k]
            gv = gold_n[k]
            if not math.isclose(pv, gv, rel_tol=rel_tol, abs_tol=abs_tol):
                diff = pv - gv
                print(
                    "[scov2_md_annotate] Value mismatch at residue "
                    f"{k}: pred={pv:.8g} gold={gv:.8g} diff={diff:.8g} "
                    f"(abs_tol={abs_tol}, rel_tol={rel_tol})"
                )
                return False

        return True

    except Exception as e:
        print(f"[scov2_md_annotate] Error processing JSON files: {e}")
        return False
