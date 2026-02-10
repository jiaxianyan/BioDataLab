import json

def m2or_annotate(d):
    try:
        # load predicted result
        with open(f'{d}/m2or_annotate.json', 'r') as f:
            pred_result = json.load(f)

        # load ground truth
        with open('benchmark/gold_results/m2or_annotate.json', 'r') as f:
            gold_result = json.load(f)

        # strict comparison (case-sensitive, key/value exact match)
        return pred_result == gold_result

    except Exception as e:
        print(f"Error processing JSON files: {e}")
        return False
