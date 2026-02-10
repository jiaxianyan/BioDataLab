import json

def amdb_retrieval(d):
    try:
        # load predicted result
        with open(f'{d}/amdb_retrieval.json', 'r') as f:
            pred_list = json.load(f)

        # load gold result
        with open('benchmark/gold_results/amdb_retrieval.json', 'r') as f:
            gold_list = json.load(f)

        # normalize: lowercase, ensure strings
        pred_list = [str(x).lower() for x in pred_list]
        gold_list = [str(x).lower() for x in gold_list]

        # order-insensitive comparison
        return sorted(pred_list) == sorted(gold_list)

    except Exception as e:
        print(f"Error processing JSON files: {e}")
        return False
