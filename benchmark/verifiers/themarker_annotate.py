import json

def themarker_annotate(d):
    try:
        # Load predicted results
        with open(f'{d}/themarker_annotate.json', 'r') as f:
            pred_labels = json.load(f)

        # Load ground truth results
        with open('benchmark/gold_results/themarker_annotate.json', 'r') as f:
            gold_labels = json.load(f)

        # Strict comparison: order-sensitive and case-sensitive
        return pred_labels == gold_labels

    except Exception as e:
        print(f"Error processing JSON files: {e}")
        return False
