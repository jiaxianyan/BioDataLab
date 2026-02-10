import csv
from collections import Counter

REQUIRED_COLS = ["chr", "start", "end", "edqtl_id", "gene"]

def _load_tsv_as_counter(path: str) -> Counter:
    """
    Load a TSV and return a multiset (Counter) of canonical 5-feature tuples:
    (chr, start, end, edqtl_id, gene)

    - Ignores row order
    - Ignores column order (as long as required columns exist)
    - Compares only the 5 key features (extra columns are ignored)
    - Normalizes:
        chr/edqtl_id/gene -> stripped + lowercased
        start/end -> int
    """
    with open(path, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"No header found in TSV: {path}")

        # Normalize header names for robust matching (e.g., "Chr" vs "chr")
        header_map = {h.strip().lower(): h for h in reader.fieldnames if h is not None}

        missing = [c for c in REQUIRED_COLS if c not in header_map]
        if missing:
            raise ValueError(
                f"Missing required columns {missing} in {path}. "
                f"Found columns: {reader.fieldnames}"
            )

        counter = Counter()
        for i, row in enumerate(reader, start=2):  # start=2 because header is line 1
            try:
                chr_val = (row[header_map["chr"]] or "").strip().lower()
                start_val = int((row[header_map["start"]] or "").strip())
                end_val = int((row[header_map["end"]] or "").strip())
                edqtl_id_val = (row[header_map["edqtl_id"]] or "").strip().lower()
                gene_val = (row[header_map["gene"]] or "").strip().lower()
            except Exception as e:
                raise ValueError(f"Bad row at line {i} in {path}: {e}. Row={row}")

            key = (chr_val, start_val, end_val, edqtl_id_val, gene_val)
            counter[key] += 1

        return counter

def gpedit_refinement(d: str) -> bool:
    """
    Verifier for GPEdit task:
    Compare predicted <d>/gpedit_refinement.tsv against gold result
    benchmark/gold_results/gpedit_refinement.tsv.

    Evaluation:
    - Must match exactly on the 5 key features: chr, start, end, edqtl_id, gene
    - Row order and column order are NOT required to match
    - Extra columns are ignored
    """
    try:
        pred_path = f"{d}/gpedit_refinement.tsv"
        gold_path = "benchmark/gold_results/gpedit_refinement.tsv"

        pred = _load_tsv_as_counter(pred_path)
        gold = _load_tsv_as_counter(gold_path)

        return pred == gold

    except FileNotFoundError as e:
        print(f"File not found: {e}")
        return False
    except Exception as e:
        print(f"Error evaluating gpedit_refinement: {e}")
        return False
