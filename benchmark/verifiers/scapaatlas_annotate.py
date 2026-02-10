import os
import csv

def scapaatlas_annotate(d: str) -> bool:
    """
    Evaluate whether the generated BED results equal the groundtruth.
    - Input:
        d: directory containing the prediction file `scapaatlas_annotate.bed`
    - Compare:
        pred: {d}/scapaatlas_annotate.bed
        gold: benchmark/gold_results/scapaatlas_annotate.bed
    - Notes:
        * order-insensitive (set comparison)
        * compare only key columns (default: first 4 columns: chrom, start, end, name)
    """
    PRED_NAME = "scapaatlas_annotate.bed"
    GOLD_PATH = os.path.join("benchmark", "gold_results", PRED_NAME)

    # Compare how many leading columns are considered "key"
    KEY_COLS = 4  # set to 3 if you only want (chrom,start,end)

    pred_path = os.path.join(d, PRED_NAME)

    try:
        if not os.path.isfile(pred_path):
            print(f"[scapaatlas_annotate] Missing prediction file: {pred_path}")
            return False
        if not os.path.isfile(GOLD_PATH):
            print(f"[scapaatlas_annotate] Missing gold file: {GOLD_PATH}")
            return False

        pred_set = _read_bed_key_set(pred_path, key_cols=KEY_COLS)
        gold_set = _read_bed_key_set(GOLD_PATH, key_cols=KEY_COLS)

        if pred_set is None or gold_set is None:
            return False

        # Exact set match
        ok = (pred_set == gold_set)
        if not ok:
            # helpful debugging
            only_in_pred = list(pred_set - gold_set)[:10]
            only_in_gold = list(gold_set - pred_set)[:10]
            print(f"[scapaatlas_annotate] Mismatch.")
            print(f"  pred rows (keyed): {len(pred_set)}")
            print(f"  gold rows (keyed): {len(gold_set)}")
            if only_in_pred:
                print(f"  examples only in pred (up to 10): {only_in_pred}")
            if only_in_gold:
                print(f"  examples only in gold (up to 10): {only_in_gold}")
        return ok

    except Exception as e:
        print(f"[scapaatlas_annotate] Error: {e}")
        return False


def _read_bed_key_set(path: str, key_cols: int = 4):
    """
    Read a BED/BED-like file into a set of tuples of key columns.
    - Ignores empty lines and comment lines starting with '#', 'track', 'browser'
    - Accepts tab/space separated or comma separated lines
    - Normalizes:
        chrom as string
        start/end as int (if present)
        name as string (if present)
    """
    key_set = set()
    skipped = 0
    parsed = 0

    def is_comment(line: str) -> bool:
        s = line.strip()
        if not s:
            return True
        lower = s.lower()
        return lower.startswith("#") or lower.startswith("track") or lower.startswith("browser")

    # Try to detect delimiter line-by-line (tabs/spaces vs commas)
    with open(path, "r", newline="") as f:
        for raw in f:
            if is_comment(raw):
                continue

            line = raw.strip()
            # Decide delimiter
            if "," in line and "\t" not in line and " " not in line:
                cols = next(csv.reader([line]))
            else:
                # split by any whitespace (handles tabs + multiple spaces)
                cols = line.split()

            if len(cols) < 3:
                skipped += 1
                continue

            # Take only key columns (but do not exceed available columns)
            use_n = min(key_cols, len(cols))
            cols = cols[:use_n]

            # Normalize required columns
            # col0: chrom
            chrom = cols[0].strip()

            # col1/2: start/end
            try:
                start = int(float(cols[1]))
                end = int(float(cols[2]))
            except Exception:
                skipped += 1
                continue

            out = [chrom, start, end]

            # col3: name (optional key col)
            if use_n >= 4:
                name = cols[3].strip()
                out.append(name)

            # If key_cols > 4 and you later want to compare more columns,
            # you can extend normalization here.
            # For now, keep any extra key columns as raw stripped strings.
            if use_n > 4:
                for x in cols[4:use_n]:
                    out.append(x.strip())

            key_set.add(tuple(out))
            parsed += 1

    if parsed == 0 and skipped > 0:
        print(f"[_read_bed_key_set] No valid rows parsed from {path} (skipped={skipped}).")
        return None

    return key_set
