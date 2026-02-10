import math

def _parse_counts_tsv(path: str) -> dict:
    """
    Parse a two-column TSV like:
      gene_id<TAB>count
    Returns {gene_id: float(count)}.
    Ignores empty lines.
    """
    counts = {}
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # Expect: ID \t count (but be tolerant to multiple spaces/tabs)
            parts = line.split("\t")
            if len(parts) < 2:
                parts = line.split()
            if len(parts) < 2:
                continue
            k = parts[0].strip()
            v = parts[1].strip()
            try:
                counts[k] = float(v)
            except ValueError:
                # If count isn't parseable, skip this row
                continue
    return counts

def _pearson_corr(x, y) -> float:
    """
    Compute Pearson correlation for two equal-length lists of numbers.
    Returns float in [-1, 1]. If either vector has zero variance, returns
    1.0 only if vectors are exactly identical, else 0.0.
    """
    n = len(x)
    if n == 0:
        return float("nan")

    mx = sum(x) / n
    my = sum(y) / n

    num = 0.0
    dx2 = 0.0
    dy2 = 0.0
    for xi, yi in zip(x, y):
        dx = xi - mx
        dy = yi - my
        num += dx * dy
        dx2 += dx * dx
        dy2 += dy * dy

    den = math.sqrt(dx2 * dy2)
    if den == 0.0:
        return 1.0 if all(xi == yi for xi, yi in zip(x, y)) else 0.0
    return num / den

def cyanoomicsdb_annotate_2(d: str) -> bool:
    """
    Verifier for cyanoomicsdb_annotate_2.

    Compares:
      Pred:  <d>/cyanoomicsdb_annotate_2.txt
      Gold:  benchmark/gold_results/cyanoomicsdb_annotate_2.txt

    Metric: Pearson correlation between predicted and gold counts over the
    UNION of feature IDs (missing treated as 0). Pass if corr > 0.98.

    Note: This includes special rows such as:
      __no_feature, __ambiguous, __too_low_aQual, __not_aligned, __alignment_not_unique
    """
    try:
        pred_path = f"{d}/cyanoomicsdb_annotate_2.txt"
        gold_path = "benchmark/gold_results/cyanoomicsdb_annotate_2.txt"

        pred = _parse_counts_tsv(pred_path)
        gold = _parse_counts_tsv(gold_path)

        # Build aligned vectors over UNION to penalize missing/extra keys.
        keys = sorted(set(pred.keys()) | set(gold.keys()))
        if len(keys) == 0:
            return False

        x = [pred.get(k, 0.0) for k in keys]
        y = [gold.get(k, 0.0) for k in keys]

        corr = _pearson_corr(x, y)
        return (not math.isnan(corr)) and (corr > 0.98)

    except Exception as e:
        print(f"Error evaluating cyanoomicsdb_annotate_2: {e}")
        return False
