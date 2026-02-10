import csv
from collections import Counter

def drmref_annotate(d):
    """
    Evaluate whether `<d>/drmref_annotate.csv` matches the ground-truth CSV
    at `benchmark/gold_results/drmref_annotate.csv`.

    Rules:
    - Case-sensitive (do NOT change gene name casing).
    - Ignore order.
    - Compare as a multiset (duplicates count).
    """
    try:
        pred_path = f"{d}/drmref_annotate.csv"
        gold_path = "benchmark/gold_results/drmref_annotate.csv"

        def read_gene_list(path):
            genes = []
            with open(path, "r", newline="", encoding="utf-8") as f:
                reader = csv.reader(f)
                rows = list(reader)

            if not rows:
                return []

            # Detect header (e.g., "gene") and decide start row
            start_idx = 0
            first_cell = rows[0][0].strip() if rows[0] else ""
            if first_cell == "gene":
                start_idx = 1

            for row in rows[start_idx:]:
                if not row:
                    continue
                g = row[0].strip()
                if g != "":
                    genes.append(g)
            return genes

        pred_genes = read_gene_list(pred_path)
        gold_genes = read_gene_list(gold_path)

        return Counter(pred_genes) == Counter(gold_genes)

    except Exception as e:
        print(f"Error processing CSV files: {e}")
        return False
