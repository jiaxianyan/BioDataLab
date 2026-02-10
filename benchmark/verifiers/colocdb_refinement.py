import csv
import math
from collections import Counter

def colocdb_refinement(d: str) -> bool:
    """
    Evaluate whether `{d}/colocdb_refinement.tsv` matches
    `benchmark/gold_results/colocdb_refinement.tsv`.

    Comparison policy:
      - TSV header must match (same column names, same set and order).
      - Row order is ignored.
      - Normalize:
          * strip whitespace
          * NA/NaN/'' -> None
          * allele columns -> lowercase
          * numeric columns -> float/int normalization (tolerant to formatting)
    """
    pred_path = f"{d}/colocdb_refinement.tsv"
    gold_path = "benchmark/gold_results/colocdb_refinement.tsv"

    # Columns that should be treated as numeric (float)
    float_cols = {
        "effect_allele_frequency", "FreqSE", "MinFreq", "MaxFreq",
        "beta", "standard_error", "p_value",
        "HetISq", "HetChiSq", "HetPVal"
    }
    # Columns that should be treated as integer
    int_cols = {"HetDf", "n", "chromosome", "base_pair_location"}

    allele_cols = {"effect_allele", "other_allele"}

    def _is_na(x: str) -> bool:
        if x is None:
            return True
        s = str(x).strip()
        return s == "" or s.upper() in {"NA", "N/A", "NAN", "NULL", "NONE"}

    def _norm_float(x: str):
        if _is_na(x):
            return None
        try:
            v = float(str(x).strip())
            if math.isnan(v):
                return None
            # rounding helps tolerate tiny representation differences
            return round(v, 12)
        except Exception:
            # if it's not parseable, keep as raw string (strict)
            return str(x).strip()

    def _norm_int(x: str):
        if _is_na(x):
            return None
        try:
            # sometimes int columns come as "1.0"
            v = float(str(x).strip())
            if math.isnan(v):
                return None
            return int(v)
        except Exception:
            return str(x).strip()

    def _norm_str(x: str):
        if _is_na(x):
            return None
        return str(x).strip()

    def _normalize_value(col: str, val: str):
        if col in allele_cols:
            v = _norm_str(val)
            return None if v is None else v.lower()
        if col in int_cols:
            return _norm_int(val)
        if col in float_cols:
            return _norm_float(val)
        # variant_id is usually case-insensitive (rsID), but keep stable:
        if col == "variant_id":
            v = _norm_str(val)
            return None if v is None else v.lower()
        return _norm_str(val)

    def _read_tsv(path: str):
        with open(path, "r", newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            if reader.fieldnames is None:
                raise ValueError(f"No header found in {path}")
            cols = [c.strip() for c in reader.fieldnames]
            rows = []
            for r in reader:
                # normalize keys in case of stray spaces
                rr = {k.strip(): v for k, v in r.items()}
                norm_row = tuple(_normalize_value(c, rr.get(c, "")) for c in cols)
                rows.append(norm_row)
            return cols, rows

    try:
        pred_cols, pred_rows = _read_tsv(pred_path)
        gold_cols, gold_rows = _read_tsv(gold_path)

        # Require identical column names and order (usually benchmark expects this)
        if pred_cols != gold_cols:
            return False

        # Compare as multisets (ignore order, but respect duplicates if any)
        return Counter(pred_rows) == Counter(gold_rows)

    except FileNotFoundError as e:
        print(f"Missing file: {e}")
        return False
    except Exception as e:
        print(f"Error processing TSV files: {e}")
        return False
