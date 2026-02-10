import json

def scqtlbase_refinement(d):
    try:
        # Load prediction
        with open(f"{d}/scqtlbase_refinement.json", "r") as f:
            pred_list = json.load(f)

        # Load ground truth
        with open("benchmark/gold_results/scqtlbase_refinement.json", "r") as f:
            gold_list = json.load(f)

        # Strict comparison:
        # - order matters
        # - case-sensitive
        # - exact match
        return pred_list == gold_list

    except Exception as e:
        print(f"Error processing JSON files: {e}")
        return False
