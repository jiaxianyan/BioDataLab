import json

def ravar_refinement_1(d):
    try:
        # load prediction
        with open(f"{d}/ravar_refinement_1.json", "r") as f:
            pred_mapping = json.load(f)

        # load ground truth
        with open("benchmark/gold_results/ravar_refinement_1.json", "r") as f:
            gold_mapping = json.load(f)

        # basic type check
        if not isinstance(pred_mapping, dict) or not isinstance(gold_mapping, dict):
            return False

        # exact dict comparison (order-independent)
        return pred_mapping == gold_mapping

    except Exception as e:
        print(f"Error processing JSON files: {e}")
        return False
