import json

def covpdb_retrieval(d):
    try:

        with open(f'{d}/cov_pdb_retrieval.json', 'r') as f:
            pred_pdb_list = json.load(f)
        # print(f'{d}/cov_pdb_retrieval.json')
        with open('benchmark/gold_results/cov_pdb_retrieval.json', 'r') as f:
            gold_pdb_list = json.load(f)

        pred_pdb_list = [pdb.lower() for pdb in pred_pdb_list]
        gold_pdb_list = [pdb.lower() for pdb in gold_pdb_list]

        return sorted(pred_pdb_list) == sorted(gold_pdb_list)
        
    except Exception as e:
        print(f"Error processing JSON files: {e}")
        return False