import json

def covid_19_extract(d):
    try:
        # load prediction
        with open(f'{d}/covid_19_extract.json', 'r') as f:
            pred_list = json.load(f)

        # load ground truth
        with open('benchmark/gold_results/covid_19_extract.json', 'r') as f:
            gold_list = json.load(f)

        # helper: normalize dict (lowercase keys and values)
        def normalize_record(record):
            return {
                str(k).lower(): str(v).lower()
                for k, v in record.items()
            }

        # normalize all records
        pred_norm = [normalize_record(r) for r in pred_list]
        gold_norm = [normalize_record(r) for r in gold_list]

        # sort lists to make comparison order-independent
        pred_norm = sorted(pred_norm, key=lambda x: json.dumps(x, sort_keys=True))
        gold_norm = sorted(gold_norm, key=lambda x: json.dumps(x, sort_keys=True))

        return pred_norm == gold_norm

    except Exception as e:
        print(f"Error processing JSON files: {e}")
        return False
