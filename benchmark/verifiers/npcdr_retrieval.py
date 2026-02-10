from .utils import load_json, evaluate_ac, AC_THRESHOLD

GT_PATH = "benchmark/gold_results/npcdr_retrieval.json"

def npcdr_retrieval(d):
    try:
        PRED_PATH = f"{d}/npcdr_retrieval.json"
        gt_lst = load_json(GT_PATH)
        pred_list = load_json(PRED_PATH)
        score = evaluate_ac(gt_lst, pred_list)
        return score >= AC_THRESHOLD
    except:
        return False