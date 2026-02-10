import os
import pandas as pd
import numpy as np

def clinicalomicsdb_annotate(d, atol=1e-6, rtol=1e-4):
    """
    Evaluate whether <d>/clinicalomicsdb_annotate.csv matches the ground-truth CSV.

    - Row/column order is NOT required to match.
    - Numeric values (p_value, auroc) are compared approximately using tolerances.
    - Requires the same set of genes and the required columns: gene, p_value, auroc.
    """
    try:
        pred_path = os.path.join(d, "clinicalomicsdb_annotate.csv")
        gold_path = os.path.join("benchmark", "gold_results", "clinicalomicsdb_annotate.csv")

        if not os.path.exists(pred_path):
            print(f"[clinicalomicsdb_annotate] Missing prediction file: {pred_path}")
            return False
        if not os.path.exists(gold_path):
            print(f"[clinicalomicsdb_annotate] Missing gold file: {gold_path}")
            return False

        pred = pd.read_csv(pred_path)
        gold = pd.read_csv(gold_path)

        # Normalize column names (strip whitespace)
        pred.columns = [c.strip() for c in pred.columns]
        gold.columns = [c.strip() for c in gold.columns]

        required = {"gene", "p_value", "auroc"}
        if not required.issubset(set(pred.columns)):
            print(f"[clinicalomicsdb_annotate] Prediction CSV missing columns: {sorted(required - set(pred.columns))}")
            return False
        if not required.issubset(set(gold.columns)):
            print(f"[clinicalomicsdb_annotate] Gold CSV missing columns: {sorted(required - set(gold.columns))}")
            return False

        # Keep only required columns
        pred = pred[list(required)].copy()
        gold = gold[list(required)].copy()

        # Normalize gene strings
        pred["gene"] = pred["gene"].astype(str).str.strip()
        gold["gene"] = gold["gene"].astype(str).str.strip()

        # Convert numeric columns robustly
        for col in ["p_value", "auroc"]:
            pred[col] = pd.to_numeric(pred[col], errors="coerce")
            gold[col] = pd.to_numeric(gold[col], errors="coerce")

        # Genes must be unique to compare cleanly
        if pred["gene"].duplicated().any():
            dups = pred.loc[pred["gene"].duplicated(), "gene"].unique()[:10]
            print(f"[clinicalomicsdb_annotate] Duplicate genes in prediction (showing up to 10): {dups}")
            return False
        if gold["gene"].duplicated().any():
            dups = gold.loc[gold["gene"].duplicated(), "gene"].unique()[:10]
            print(f"[clinicalomicsdb_annotate] Duplicate genes in gold (showing up to 10): {dups}")
            return False

        pred = pred.set_index("gene").sort_index()
        gold = gold.set_index("gene").sort_index()

        # Same gene set (row order not required)
        if set(pred.index) != set(gold.index):
            missing_in_pred = sorted(set(gold.index) - set(pred.index))[:10]
            extra_in_pred = sorted(set(pred.index) - set(gold.index))[:10]
            print("[clinicalomicsdb_annotate] Gene set mismatch.")
            if missing_in_pred:
                print(f"  Missing in prediction (up to 10): {missing_in_pred}")
            if extra_in_pred:
                print(f"  Extra in prediction (up to 10): {extra_in_pred}")
            return False

        # Align exactly by gene
        pred = pred.loc[gold.index]

        # Compare numeric columns with tolerance; treat NaN vs NaN as equal
        for col in ["p_value", "auroc"]:
            a = pred[col].to_numpy(dtype=float)
            b = gold[col].to_numpy(dtype=float)

            both_nan = np.isnan(a) & np.isnan(b)
            close = np.isclose(a, b, atol=atol, rtol=rtol, equal_nan=True)

            # equal_nan=True already treats NaN==NaN as True, but keep explicit for clarity
            ok = close | both_nan
            if not np.all(ok):
                bad_idx = np.where(~ok)[0][:10]
                bad_genes = gold.index[bad_idx].tolist()
                print(f"[clinicalomicsdb_annotate] Value mismatch in column '{col}' (showing up to 10 genes): {bad_genes}")
                # Print a small diff sample
                for g in bad_genes[:5]:
                    print(f"  {g}: pred={pred.loc[g, col]} gold={gold.loc[g, col]}")
                return False

        return True

    except Exception as e:
        print(f"[clinicalomicsdb_annotate] Error: {e}")
        return False
