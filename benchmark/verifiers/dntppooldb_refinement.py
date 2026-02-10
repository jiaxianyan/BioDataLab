import os
from collections import Counter

import pandas as pd


def dntppooldb_refinement(d: str) -> bool:
    """
    Verifier for dntppooldb_refinement.

    Compares:
      - prediction: {d}/dntppooldb_refinement.csv
      - ground truth: benchmark/gold_results/dntppooldb_refinement.csv

    Rules:
      - Row order does NOT matter
      - Column order does NOT matter
      - Numeric values can be approximately equal (handled via rounding)
    """
    try:
        pred_path = os.path.join(d, "dntppooldb_refinement.csv")
        gold_path = os.path.join("benchmark", "gold_results", "dntppooldb_refinement.csv")

        if not os.path.exists(pred_path):
            print(f"Missing prediction file: {pred_path}")
            return False
        if not os.path.exists(gold_path):
            print(f"Missing gold file: {gold_path}")
            return False

        pred_df = pd.read_csv(pred_path)
        gold_df = pd.read_csv(gold_path)

        # Normalize column names (strip spaces)
        pred_df.columns = [str(c).strip() for c in pred_df.columns]
        gold_df.columns = [str(c).strip() for c in gold_df.columns]

        pred_cols = set(pred_df.columns)
        gold_cols = set(gold_df.columns)

        if pred_cols != gold_cols:
            print(f"Column mismatch.\nPred: {sorted(pred_cols)}\nGold: {sorted(gold_cols)}")
            return False

        # Reorder columns deterministically (order doesn't matter)
        col_order = sorted(gold_cols)
        pred_df = pred_df[col_order].copy()
        gold_df = gold_df[col_order].copy()

        # Normalize cell values:
        # - Strip whitespace for strings
        # - Convert numeric-like columns to numeric and round to allow approximate comparison
        def normalize_df(df: pd.DataFrame) -> pd.DataFrame:
            out = df.copy()

            # Strip whitespace for object columns
            for c in out.columns:
                if out[c].dtype == object:
                    out[c] = out[c].apply(lambda x: x.strip() if isinstance(x, str) else x)

            # Detect numeric columns by coercion ratio
            for c in out.columns:
                coerced = pd.to_numeric(out[c], errors="coerce")
                non_na = out[c].notna().sum()
                num_na = coerced.notna().sum()

                # If most non-null values are numeric-like, treat as numeric
                if non_na == 0:
                    continue
                if num_na / non_na >= 0.9:
                    out[c] = pd.to_numeric(out[c], errors="coerce").round(6)

            # Make missing values comparable
            return out.where(out.notna(), other="__NA__")

        pred_n = normalize_df(pred_df)
        gold_n = normalize_df(gold_df)

        # Compare as multisets of rows (row order ignored)
        pred_rows = Counter(tuple(row) for row in pred_n.itertuples(index=False, name=None))
        gold_rows = Counter(tuple(row) for row in gold_n.itertuples(index=False, name=None))

        if pred_rows != gold_rows:
            print("Row content mismatch (order ignored; numeric approx via rounding).")
            return False

        return True

    except Exception as e:
        print(f"Error processing CSV files: {e}")
        return False
