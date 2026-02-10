import json

def vimic_extract(d):
    try:
        # Load predicted result
        with open(f'{d}/vimic_extract.json', 'r') as f:
            pred_mutations = json.load(f)

        # Load ground truth
        with open('benchmark/gold_results/vimic_extract.json', 'r') as f:
            gold_mutations = json.load(f)

        # Ensure both are lists
        if not isinstance(pred_mutations, list) or not isinstance(gold_mutations, list):
            return False

        # Order-insensitive, case-sensitive comparison
        return sorted(pred_mutations) == sorted(gold_mutations)

    except Exception as e:
        print(f"Error processing JSON files: {e}")
        return False
