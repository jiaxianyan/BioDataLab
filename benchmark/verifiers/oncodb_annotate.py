import os
import pandas as pd
import numpy as np

def oncodb_annotate(d: str) -> bool:
    """
    Evaluate whether generated oncodb_annotate.csv matches the ground truth CSV.

    - Input:
        d: the run/output directory root (where <r></r>/oncodb_annotate.csv is saved)

    - Logic:
        * Compare key columns only (covariate + major numeric stats).
        * Row/column order not required.
        * Numeric values compared approximately with tolerances.
    """
    try:
        pred_path = os.path.join(d, "oncodb_annotate.csv")
        gold_path = os.path.join("benchmark", "gold_results", "oncodb_annotate.csv")

        if not os.path.exists(pred_path):
            print(f"Missing prediction file: {pred_path}")
            return False
        if not os.path.exists(gold_path):
            print(f"Missing gold file: {gold_path}")
            return False

        pred = pd.read_csv(pred_path)
        gold = pd.read_csv(gold_path)

        # Normalize column names (strip spaces)
        pred.columns = [str(c).strip() for c in pred.columns]
        gold.columns = [str(c).strip() for c in gold.columns]

        if "covariate" not in pred.columns or "covariate" not in gold.columns:
            print("Missing required 'covariate' column.")
            return False

        # Normalize covariate keys
        pred["covariate"] = pred["covariate"].astype(str).str.strip().str.lower()
        gold["covariate"] = gold["covariate"].astype(str).str.strip().str.lower()

        # Keep only essential columns (if present)
        key_cols = [
            "coef",
            "exp(coef)",
            "se(coef)",
            "z",
            "p",
            "-log2(p)",
            "coef lower 95%",
            "coef upper 95%",
            "exp(coef) lower 95%",
            "exp(coef) upper 95%",
        ]
        pred_keep = ["covariate"] + [c for c in key_cols if c in pred.columns]
        gold_keep = ["covariate"] + [c for c in key_cols if c in gold.columns]

        # Compare on intersection of available numeric columns
        common_numeric = [c for c in key_cols if (c in pred.columns and c in gold.columns)]
        if len(common_numeric) == 0:
            # Still allow a "structure-only" check: covariate set must match
            pred_cov = set(pred["covariate"].dropna().tolist())
            gold_cov = set(gold["covariate"].dropna().tolist())
            return pred_cov == gold_cov

        pred_s = pred[pred_keep].copy()
        gold_s = gold[gold_keep].copy()

        # Deduplicate by covariate (keep first) to avoid merge explosions
        pred_s = pred_s.drop_duplicates(subset=["covariate"], keep="first")
        gold_s = gold_s.drop_duplicates(subset=["covariate"], keep="first")

        # Covariate coverage must match (strict)
        pred_cov = set(pred_s["covariate"].dropna().tolist())
        gold_cov = set(gold_s["covariate"].dropna().tolist())
        if pred_cov != gold_cov:
            # If you want to be slightly lenient, you can relax this,
            # but benchmark tasks usually expect same covariates.
            missing_in_pred = sorted(list(gold_cov - pred_cov))
            extra_in_pred = sorted(list(pred_cov - gold_cov))
            print(f"Covariate mismatch. Missing in pred: {missing_in_pred}, Extra in pred: {extra_in_pred}")
            return False

        # Merge aligned by covariate (row order irrelevant)
        m = gold_s[["covariate"] + common_numeric].merge(
            pred_s[["covariate"] + common_numeric],
            on="covariate",
            how="inner",
            suffixes=("_gold", "_pred"),
        )

        # Convert numeric columns robustly
        for c in common_numeric:
            m[f"{c}_gold"] = pd.to_numeric(m[f"{c}_gold"], errors="coerce")
            m[f"{c}_pred"] = pd.to_numeric(m[f"{c}_pred"], errors="coerce")

        # Tolerances: allow small numeric drift (float formatting / library differences)
        # - absolute tol: small values like p might differ a bit
        # - relative tol: larger coefficients can differ proportionally
        abs_tol = 1e-3
        rel_tol = 2e-2  # 2%

        # Evaluate column-by-column; allow a small fraction of failures overall
        total_checks = 0
        total_fail = 0

        for c in common_numeric:
            g = m[f"{c}_gold"].to_numpy(dtype=float)
            p = m[f"{c}_pred"].to_numpy(dtype=float)

            # If both are NaN for a row, treat as pass; if one NaN other not, fail
            both_nan = np.isnan(g) & np.isnan(p)
            one_nan = np.isnan(g) ^ np.isnan(p)
            if np.any(one_nan):
                total_checks += int(np.sum(one_nan))
                total_fail += int(np.sum(one_nan))

            # For non-nan pairs, compare approximately
            mask = ~np.isnan(g) & ~np.isnan(p)
            if np.any(mask):
                total_checks += int(np.sum(mask))
                ok = np.isclose(p[mask], g[mask], rtol=rel_tol, atol=abs_tol)

                # Special handling for p-values (often tiny; allow a bit more relative slack)
                if c == "p":
                    ok = np.isclose(p[mask], g[mask], rtol=5e-2, atol=1e-4)

                total_fail += int(np.sum(~ok))

        if total_checks == 0:
            # Nothing numeric was actually comparable; fall back to covariate set check already done
            return True

        fail_rate = total_fail / total_checks

        # Accept if most key numbers match
        # You can tune this threshold depending on how strict you want to be.
        return fail_rate <= 0.05  # <= 5% mismatches

    except Exception as e:
        print(f"Error evaluating oncodb_annotate: {e}")
        return False
