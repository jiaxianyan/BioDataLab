import json

def diana_mited_retrieval(d):
    """
    Evaluate whether the generated accession list matches the ground truth.

    Args:
        d (str): The root directory where the result json is saved.

    Returns:
        bool: True if prediction equals ground truth (case-sensitive, order-sensitive),
              False otherwise.
    """
    try:
        # load predicted result
        with open(f'{d}/diana_mited_retrieval.json', 'r') as f:
            pred_list = json.load(f)

        # load gold result
        with open('benchmark/gold_results/diana_mited_retrieval.json', 'r') as f:
            gold_list = json.load(f)

        # directly compare lists
        return pred_list == gold_list

    except Exception as e:
        print(f"Error processing JSON files: {e}")
        return False
