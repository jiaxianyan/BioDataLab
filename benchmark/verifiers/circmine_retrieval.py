from .utils import load_json, evaluate_ac, AC_THRESHOLD

GT_PATH = "benchmark/gold_results/circmine_retrieval.json"

def circmine_retrieval(d):
    try:
        PRED_PATH = f"{d}/circmine_retrieval.json"
        gt_lst = load_json(GT_PATH)
        pred_list = load_json(PRED_PATH)
        score = evaluate_ac(gt_lst, pred_list,"Series")
        return score >= AC_THRESHOLD
    except:
        return False