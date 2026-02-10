import os
import re
import math
import pandas as pd

def covpdb_integration(d,
             pred_name="covpdb_integration_complex.csv",
             gold_path="benchmark/gold_results/covpdb_integration_complex.csv",
             rtol=1e-3,
             atol=1e-2,
             affinity_rtol=1e-2,
             affinity_atol=1e-2):
    """
    Evaluate whether predicted csv equals groundtruth csv with approximate equality for numeric columns.

    Args:
        d (str): directory containing prediction csv (d/covpdb_integration_complex.csv)
        pred_name (str): prediction filename
        gold_path (str): ground truth csv path
        rtol/atol: tolerance for Resolution and SASA numeric comparisons
        affinity_rtol/affinity_atol: tolerance for Affinity numeric comparisons (after parsing)

    Returns:
        bool
    """
    try:
        pred_path = os.path.join(d, pred_name)
        if not os.path.exists(pred_path):
            print(f"[EVAL] Pred file not found: {pred_path}")
            return False
        if not os.path.exists(gold_path):
            print(f"[EVAL] Gold file not found: {gold_path}")
            return False

        pred = pd.read_csv(pred_path)
        gold = pd.read_csv(gold_path)

        # ---- basic schema checks ----
        required_cols = ["PDB_ID", "Method", "Resolution", "Affinity", "SASA"]
        for c in required_cols:
            if c not in pred.columns:
                print(f"[EVAL] Missing column in pred: {c}")
                return False
            if c not in gold.columns:
                print(f"[EVAL] Missing column in gold: {c}")
                return False

        # ---- normalize helpers ----
        def norm_text(x):
            if pd.isna(x):
                return ""
            return str(x).strip().lower()

        def norm_method(x):
            s = norm_text(x)
            # light normalization: remove extra spaces/punct variations
            s = re.sub(r"\s+", " ", s)
            s = s.replace("_", " ").replace("-", " ")
            s = re.sub(r"\s+", " ", s).strip()
            return s

        def to_float(x):
            if pd.isna(x):
                return None
            if isinstance(x, (int, float)) and not (isinstance(x, float) and math.isnan(x)):
                return float(x)
            s = str(x).strip()
            # extract first float-like token
            m = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", s)
            if not m:
                return None
            try:
                return float(m.group(0))
            except:
                return None

        def is_close(a, b, rtol_, atol_):
            if a is None or b is None:
                return False
            return math.isclose(a, b, rel_tol=rtol_, abs_tol=atol_)

        def norm_affinity_string(s):
            # fallback: ignore case and repeated spaces
            s = norm_text(s)
            s = re.sub(r"\s+", " ", s)
            return s

        def parse_affinity(val):
            """
            Try parse affinity strings like:
              "IC50 50 nM", "Kd=1.2 uM", "Ki: 10 nM", "EC50 0.5 mM"
            Returns (kind, number, unit) or None.
            """
            if pd.isna(val):
                return None
            s = str(val).strip()
            s = re.sub(r"\s+", " ", s)

            # kind
            kind_m = re.search(r"\b(ic50|ki|kd|ec50|kon|koff)\b", s, flags=re.IGNORECASE)
            kind = kind_m.group(1).lower() if kind_m else None

            # number
            num_m = re.search(r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", s)
            num = float(num_m.group(1)) if num_m else None

            # unit (nm, um, mm, m, pm; allow "nM" "uM" "µM" "mM")
            unit_m = re.search(r"\b(pM|nM|uM|µM|mM|M)\b", s)
            unit = unit_m.group(1) if unit_m else None
            if unit == "µM":
                unit = "uM"

            if kind is None and num is None and unit is None:
                return None
            if num is None:
                return None
            return (kind if kind else "", num, unit if unit else "")

        # ---- align rows by PDB_ID (case-insensitive) ----
        pred["_pdb_key"] = pred["PDB_ID"].apply(norm_text)
        gold["_pdb_key"] = gold["PDB_ID"].apply(norm_text)

        pred_keys = set(pred["_pdb_key"].tolist())
        gold_keys = set(gold["_pdb_key"].tolist())

        if pred_keys != gold_keys:
            missing = sorted(list(gold_keys - pred_keys))
            extra = sorted(list(pred_keys - gold_keys))
            if missing:
                print(f"[EVAL] Missing PDB_IDs in pred: {missing[:20]}{' ...' if len(missing)>20 else ''}")
            if extra:
                print(f"[EVAL] Extra PDB_IDs in pred: {extra[:20]}{' ...' if len(extra)>20 else ''}")
            return False

        # build lookup
        pred_map = pred.set_index("_pdb_key")
        gold_map = gold.set_index("_pdb_key")

        # ---- compare row by row ----
        for k in sorted(list(gold_keys)):
            pr = pred_map.loc[k]
            gr = gold_map.loc[k]

            # if duplicate rows for same pdb_id, fail (or you can choose to handle)
            if isinstance(pr, pd.DataFrame) or isinstance(gr, pd.DataFrame):
                print(f"[EVAL] Duplicate PDB_ID rows detected for {k}.")
                return False

            # Method (strict after normalization)
            if norm_method(pr["Method"]) != norm_method(gr["Method"]):
                print(f"[EVAL] Method mismatch for {k}: pred='{pr['Method']}' gold='{gr['Method']}'")
                return False

            # Resolution (approx)
            pr_res = to_float(pr["Resolution"])
            gr_res = to_float(gr["Resolution"])
            if not is_close(pr_res, gr_res, rtol, atol):
                print(f"[EVAL] Resolution mismatch for {k}: pred={pr['Resolution']} gold={gr['Resolution']}")
                return False

            # SASA (approx)
            pr_sasa = to_float(pr["SASA"])
            gr_sasa = to_float(gr["SASA"])
            if not is_close(pr_sasa, gr_sasa, rtol, atol):
                print(f"[EVAL] SASA mismatch for {k}: pred={pr['SASA']} gold={gr['SASA']}")
                return False

            # Affinity (parse + approx; fallback to normalized string)
            pr_aff_p = parse_affinity(pr["Affinity"])
            gr_aff_p = parse_affinity(gr["Affinity"])

            if pr_aff_p and gr_aff_p:
                pr_kind, pr_num, pr_unit = pr_aff_p
                gr_kind, gr_num, gr_unit = gr_aff_p
                if pr_kind != gr_kind:
                    print(f"[EVAL] Affinity kind mismatch for {k}: pred='{pr_kind}' gold='{gr_kind}'")
                    return False
                if pr_unit.lower() != gr_unit.lower():
                    print(f"[EVAL] Affinity unit mismatch for {k}: pred='{pr_unit}' gold='{gr_unit}'")
                    return False
                if not is_close(pr_num, gr_num, affinity_rtol, affinity_atol):
                    print(f"[EVAL] Affinity value mismatch for {k}: pred={pr['Affinity']} gold={gr['Affinity']}")
                    return False
            else:
                # fallback: string compare after normalization
                if norm_affinity_string(pr["Affinity"]) != norm_affinity_string(gr["Affinity"]):
                    print(f"[EVAL] Affinity mismatch for {k}: pred='{pr['Affinity']}' gold='{gr['Affinity']}'")
                    return False

        return True

    except Exception as e:
        print(f"[EVAL] Error evaluating CSVs: {e}")
        return False
