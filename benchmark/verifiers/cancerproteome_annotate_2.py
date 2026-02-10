import os

def cancerproteome_annotate_2(d):
    """
    Evaluate whether the generated cancerproteome_annotate_2.txt matches the groundtruth.
    - Input: d (str) directory path that contains the predicted result file.
    - Prediction file: <d>/cancerproteome_annotate_2.txt
    - Groundtruth file: benchmark/gold_results/cancerproteome_annotate_2.txt
    Rules:
    - Each non-header line is a protein pair: ProteinA,ProteinB
    - Pair order in file is irrelevant
    - Within a pair, (A,B) is considered the same as (B,A)
    - Ignore empty lines and tolerate presence/absence of a header row.
    """
    pred_path = os.path.join(d, "cancerproteome_annotate_2.txt")
    gold_path = os.path.join("benchmark", "gold_results", "cancerproteome_annotate_2.txt")

    def _read_pairs(path):
        if not os.path.exists(path):
            raise FileNotFoundError(f"Missing file: {path}")

        pairs = set()
        with open(path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                # Tolerate a header row like: ProteinA,ProteinB (case-insensitive, spaces allowed)
                header_norm = "".join(line.lower().split())
                if header_norm in ("proteina,proteinb", "protein_a,protein_b"):
                    continue

                parts = [p.strip() for p in line.split(",")]
                if len(parts) != 2:
                    # If there are malformed lines, treat as evaluation failure
                    raise ValueError(f"Malformed line (expected 2 comma-separated fields): {line}")

                a, b = parts[0], parts[1]
                if not a or not b:
                    raise ValueError(f"Empty protein name in line: {line}")

                # Canonicalize pair so (A,B) == (B,A); keep original case but canonical sort
                # Also normalize internal whitespace
                a_norm = " ".join(a.split())
                b_norm = " ".join(b.split())
                pair = tuple(sorted((a_norm, b_norm), key=lambda x: x.lower()))
                pairs.add(pair)

        return pairs

    try:
        pred_pairs = _read_pairs(pred_path)
        gold_pairs = _read_pairs(gold_path)

        return pred_pairs == gold_pairs

    except Exception as e:
        print(f"Error evaluating cancerproteome_annotate_2: {e}")
        return False
