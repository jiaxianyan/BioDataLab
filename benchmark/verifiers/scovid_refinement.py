import os
import csv
from collections import Counter

def scovid_refinement(d: str) -> bool:
    """
    Evaluate whether the predicted qc_passed_barcodes.csv matches the gold result.

    Pred file:  {d}/qc_passed_barcodes.csv
    Gold file:  benchmark/gold_results/scovid_refinement.csv

    Comparison:
      - Ignore row order
      - Compare as multiset (Counter) to catch duplicates/missing/extra rows
      - Strip whitespace, drop empty lines
      - Treat each row as a single barcode string (first column if multiple columns)
    """
    try:
        pred_path = os.path.join(d, "qc_passed_barcodes.csv")
        gold_path = 'benchmark/gold_results/scovid_refinement.csv'

        if not os.path.exists(pred_path):
            print(f"[scov2_md_annotate] Missing prediction file: {pred_path}")
            return False
        if not os.path.exists(gold_path):
            print(f"[scov2_md_annotate] Missing gold file: {gold_path}")
            return False

        def load_barcodes(path: str):
            items = []
            with open(path, "r", newline="") as f:
                reader = csv.reader(f)
                for row in reader:
                    if not row:
                        continue
                    # task要求 single column；但这里更鲁棒：取第一列即可
                    s = row[0].strip()
                    if s:
                        items.append(s)
            return items

        pred = load_barcodes(pred_path)
        gold = load_barcodes(gold_path)

        # 如果你希望“大小写不敏感”，可以改成 s.lower()；但一般 barcode/GSM 大小写应严格一致
        pred_counter = Counter(pred)
        gold_counter = Counter(gold)

        ok = (pred_counter == gold_counter)

        if not ok:
            pred_set, gold_set = set(pred_counter.keys()), set(gold_counter.keys())
            missing = sorted(list(gold_set - pred_set))[:10]
            extra = sorted(list(pred_set - gold_set))[:10]
            print(f"[scov2_md_annotate] Mismatch.")
            print(f"  pred rows: {len(pred)} (unique {len(pred_set)})")
            print(f"  gold rows: {len(gold)} (unique {len(gold_set)})")
            if missing:
                print(f"  missing (first 10): {missing}")
            if extra:
                print(f"  extra   (first 10): {extra}")

            # 检查重复行（如果有）
            pred_dups = [k for k, v in pred_counter.items() if v > 1]
            if pred_dups:
                print(f"  pred has duplicates (first 10): {pred_dups[:10]}")

        return ok

    except Exception as e:
        print(f"[scov2_md_annotate] Error: {e}")
        return False
