def fusionneoantigen_annotate_2(d: str) -> bool:
    """
    Evaluate whether the predicted result matches the ground truth for FusionNeoAntigen docking task.

    - Pred file:   f"{d}/fusionneoantigen_annotate_2.txt"
    - Gold file:   "benchmark/gold_results/fusionneoantigen_annotate_2.txt"
    - Logic: compare the (non-empty) SMILES line(s). Typically it's a single line; we compare
      the cleaned non-empty lines exactly.
    """
    try:
        pred_path = f"{d}/fusionneoantigen_annotate_2.txt"
        gold_path = "benchmark/gold_results/fusionneoantigen_annotate_2.txt"

        def read_nonempty_lines(path: str):
            with open(path, "r", encoding="utf-8") as f:
                return [ln.strip() for ln in f.read().splitlines() if ln.strip()]

        pred_lines = read_nonempty_lines(pred_path)
        gold_lines = read_nonempty_lines(gold_path)

        # Task spec implies exactly one SMILES line; still handle extra empty lines safely.
        # If multiple non-empty lines exist, require exact match in order/content.
        return pred_lines == gold_lines

    except Exception as e:
        print(f"Error evaluating fusionneoantigen_annotate_2: {e}")
        return False
