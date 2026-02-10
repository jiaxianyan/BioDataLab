def diana_mited_refinement(d, threshold=0.5):
    """
    Evaluate whether the predicted miRNA assignment ratio matches the ground truth
    within an acceptable error threshold.

    Args:
        d (str): directory containing the predicted result file
        threshold (float): allowed absolute difference in percentage points

    Returns:
        bool: True if prediction is within threshold, False otherwise
    """
    try:
        # Predicted result
        pred_path = f"{d}/diana_mited_refinement.txt"
        with open(pred_path, "r") as f:
            pred_raw = f.read().strip()

        # Gold result
        gold_path = "benchmark/gold_results/diana_mited_refinement.txt"
        with open(gold_path, "r") as f:
            gold_raw = f.read().strip()

        # Parse percentage values
        pred_value = float(pred_raw.replace("%", "").strip())
        gold_value = float(gold_raw.replace("%", "").strip())

        # Compare with threshold
        return abs(pred_value - gold_value) <= threshold

    except Exception as e:
        print(f"Error evaluating diana_mited_refinement: {e}")
        return False
