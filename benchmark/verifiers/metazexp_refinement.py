import json

def metazexp_refinement(d):
    try:
        # Load predicted result
        with open(f'{d}/metazexp_refinement.json', 'r') as f:
            pred_list = json.load(f)

        # Load gold result
        with open('benchmark/gold_results/metazexp_refinement.json', 'r') as f:
            gold_list = json.load(f)

        # Normalize: ensure list of strings, strip spaces
        pred_list = [str(x).strip() for x in pred_list]
        gold_list = [str(x).strip() for x in gold_list]

        # Compare ignoring order
        return sorted(pred_list) == sorted(gold_list)

    except Exception as e:
        print(f"Error processing metazexp_refinement evaluation: {e}")
        return False
