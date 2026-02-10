import json

def pcmdb_extract_2(d):
    try:
        # load prediction
        with open(f'{d}/pcmdb_extract_2.json', 'r') as f:
            pred_list = json.load(f)

        # load ground truth
        with open('benchmark/gold_results/pcmdb_extract_2.json', 'r') as f:
            gold_list = json.load(f)

        # normalize: lowercase, strip whitespace
        pred_list = [p.strip().lower() for p in pred_list]
        gold_list = [g.strip().lower() for g in gold_list]

        # order-insensitive comparison
        return sorted(pred_list) == sorted(gold_list)

    except Exception as e:
        print(f"Error processing JSON files: {e}")
        return False
