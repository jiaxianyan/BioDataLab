from .utils import load_json, evaluate_lst, CVS_THRESHOLD

GT_PATH = "benchmark/gold_results/asmdb_retrieval.json"

def asmdb_retrieval(d):
    try:
        PRED_PATH = f"{d}/asmdb_retrieval.json"
        gt_lst = load_json(GT_PATH)
        pred_list = load_json(PRED_PATH)
        score = evaluate_lst(gt_lst, pred_list)
        return score >= CVS_THRESHOLD
    except:
        return False