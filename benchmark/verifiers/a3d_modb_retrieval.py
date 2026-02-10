import json

def a3d_modb_retrieval(d):
    try:
        # load predicted result
        with open(f'{d}/a3d_modb_retrieval.json', 'r') as f:
            pred_list = json.load(f)

        # load gold result
        with open('benchmark/gold_results/a3d_modb_retrieval.json', 'r') as f:
            gold_list = json.load(f)

        # normalize UniProt accessions:
        # strip whitespace and convert to uppercase
        pred_list = [acc.strip().upper() for acc in pred_list]
        gold_list = [acc.strip().upper() for acc in gold_list]

        # compare ignoring order
        return sorted(pred_list) == sorted(gold_list)

    except Exception as e:
        print(f"Error processing JSON files: {e}")
        return False
