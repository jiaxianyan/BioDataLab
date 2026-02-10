import json

def tf_marker_annotate(d):
    try:
        # load prediction
        with open(f'{d}/tf_marker_annotate.json', 'r') as f:
            pred_list = json.load(f)

        # load ground truth
        with open('benchmark/gold_results/tf_marker_annotate.json', 'r') as f:
            gold_list = json.load(f)

        # ensure both are lists
        if not isinstance(pred_list, list) or not isinstance(gold_list, list):
            return False

        # case-insensitive comparison, order-sensitive
        pred_list = [str(x).lower() for x in pred_list]
        gold_list = [str(x).lower() for x in gold_list]

        return pred_list == gold_list

    except Exception as e:
        print(f"Error processing tf_marker_annotate.json: {e}")
        return False
