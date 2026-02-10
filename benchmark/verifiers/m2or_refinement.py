def m2or_refinement(d):
    """
    Evaluate whether the generated InChIKey equals the ground truth InChIKey.
    
    Args:
        d (str): directory path that contains the model output file
    
    Returns:
        bool: True if prediction matches ground truth, False otherwise
    """
    try:
        # path to predicted result
        pred_path = f"{d}/m2or_refinement.txt"
        # path to ground truth result
        gold_path = "benchmark/gold_results/m2or_refinement.txt"

        with open(pred_path, "r") as f:
            pred_inchikey = f.read().strip()

        with open(gold_path, "r") as f:
            gold_inchikey = f.read().strip()

        # Normalize (InChIKey is case-insensitive but usually uppercase)
        pred_inchikey = pred_inchikey.upper()
        gold_inchikey = gold_inchikey.upper()

        return pred_inchikey == gold_inchikey

    except Exception as e:
        print(f"Error during m2or_refinement evaluation: {e}")
        return False
