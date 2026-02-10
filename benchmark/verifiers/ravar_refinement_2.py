import json

def ravar_refinement_2(d):
    try:
        # Load prediction
        with open(f'{d}/ravar_refinement_2.json', 'r') as f:
            pred_list = json.load(f)

        # Load ground truth
        with open('benchmark/gold_results/ravar_refinement_2.json', 'r') as f:
            gold_list = json.load(f)

        # Ensure both are lists
        if not isinstance(pred_list, list) or not isinstance(gold_list, list):
            return False

        # Order-sensitive & case-sensitive comparison
        return pred_list == gold_list

    except Exception as e:
        print(f"Error processing RAVAR refinement results: {e}")
        return False
