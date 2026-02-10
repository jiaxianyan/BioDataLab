import json

def ddinter_annotate_1(d):
    try:
        # load predicted result
        with open(f'{d}/ddinter_annotate_1.json', 'r') as f:
            pred_list = json.load(f)

        # load ground truth
        with open('benchmark/gold_results/ddinter_annotate_1.json', 'r') as f:
            gold_list = json.load(f)

        # case-insensitive comparison, order-sensitive
        pred_list = [x.lower() for x in pred_list]
        gold_list = [x.lower() for x in gold_list]

        return pred_list == gold_list

    except Exception as e:
        print(f"Error processing JSON files: {e}")
        return False
