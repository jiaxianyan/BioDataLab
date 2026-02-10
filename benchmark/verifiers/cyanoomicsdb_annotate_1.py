import csv
import os

EXPECTED_COLS = [
    "Locus_Tag",
    "Symbol",
    "Old_Locus_Tag",
    "Location",
    "Protein_ID",
    "Product",
    "Pfam",
]

def _norm_cell(x: str) -> str:
    if x is None:
        return "-"
    x = str(x).strip()
    return x if x != "" else "-"

def _norm_pfam(x: str) -> str:
    x = _norm_cell(x)
    if x == "-" or x.lower() == "na":
        return "-"
    parts = [p.strip() for p in x.split(";") if p.strip()]
    if not parts:
        return "-"
    # sort + deduplicate (stable)
    parts = sorted(set(parts))
    return ";".join(parts)

def _load_tsv(path: str):
    with open(path, "r", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        rows = list(reader)

    if not rows:
        raise ValueError(f"Empty TSV: {path}")

    header = [h.strip() for h in rows[0]]
    header_map = {h: i for i, h in enumerate(header)}

    # Require all expected columns to exist
    for col in EXPECTED_COLS:
        if col not in header_map:
            raise ValueError(f"Missing column '{col}' in {path}. Found: {header}")

    data = {}
    for r in rows[1:]:
        if not r or all(str(c).strip() == "" for c in r):
            continue

        def get(col):
            i = header_map[col]
            return r[i] if i < len(r) else "-"

        locus = _norm_cell(get("Locus_Tag"))
        if locus == "-":
            # invalid row
            continue

        rec = {
            "Locus_Tag": locus,
            "Symbol": _norm_cell(get("Symbol")),
            "Old_Locus_Tag": _norm_cell(get("Old_Locus_Tag")),
            "Location": _norm_cell(get("Location")),
            "Protein_ID": _norm_cell(get("Protein_ID")),
            "Product": _norm_cell(get("Product")),
            "Pfam": _norm_pfam(get("Pfam")),
        }

        # If duplicate locus tags appear, keep the first and treat mismatch as failure later
        if locus in data and data[locus] != rec:
            # mark a special collision
            data[locus] = ("__DUPLICATE_CONFLICT__", data[locus], rec)
        else:
            data[locus] = rec

    return data

def cyanoomicsdb_annotate_1(d: str) -> bool:
    """
    Evaluate whether the generated TSV equals the groundtruth TSV.

    Pred path:  {d}/cyanoomicsdb_annotate_1.tsv
    Gold path:  benchmark/gold_results/cyanoomicsdb_annotate_1.tsv
    """
    try:
        pred_path = os.path.join(d, "cyanoomicsdb_annotate_1.tsv")
        gold_path = "benchmark/gold_results/cyanoomicsdb_annotate_1.tsv"

        if not os.path.exists(pred_path):
            print(f"Missing prediction file: {pred_path}")
            return False
        if not os.path.exists(gold_path):
            print(f"Missing gold file: {gold_path}")
            return False

        pred = _load_tsv(pred_path)
        gold = _load_tsv(gold_path)

        # Fail on duplicate conflicts
        for k, v in pred.items():
            if isinstance(v, tuple) and v and v[0] == "__DUPLICATE_CONFLICT__":
                print(f"Duplicate locus tag with conflicting rows in prediction: {k}")
                return False
        for k, v in gold.items():
            if isinstance(v, tuple) and v and v[0] == "__DUPLICATE_CONFLICT__":
                print(f"Duplicate locus tag with conflicting rows in gold: {k}")
                return False

        # Compare keys
        if set(pred.keys()) != set(gold.keys()):
            missing = sorted(set(gold.keys()) - set(pred.keys()))
            extra = sorted(set(pred.keys()) - set(gold.keys()))
            if missing:
                print(f"Missing locus tags (in pred): {missing[:10]}{' ...' if len(missing) > 10 else ''}")
            if extra:
                print(f"Extra locus tags (in pred): {extra[:10]}{' ...' if len(extra) > 10 else ''}")
            return False

        # Compare row-by-row (Pfam normalized; other fields stripped)
        for locus in gold.keys():
            if pred[locus] != gold[locus]:
                print(f"Mismatch at locus tag: {locus}")
                print("Pred:", pred[locus])
                print("Gold:", gold[locus])
                return False

        return True

    except Exception as e:
        print(f"Error processing TSV files: {e}")
        return False
