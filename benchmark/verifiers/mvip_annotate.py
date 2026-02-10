import json

def mvip_annotate(d):
    try:
        # load prediction
        with open(f'{d}/mvip_annotate.tsv', 'r') as f:
            pred_idx_list = json.load(f)

        # load gold
        with open('benchmark/gold_results/mvip_annotate.json', 'r') as f:
            gold_idx_list = json.load(f)

        # normalize: convert all indices to int
        pred_idx_list = [int(i) for i in pred_idx_list]
        gold_idx_list = [int(i) for i in gold_idx_list]

        # compare ignoring order
        return sorted(pred_idx_list) == sorted(gold_idx_list)

    except Exception as e:
        print(f"Error processing mvip_annotate results: {e}")
        return False
