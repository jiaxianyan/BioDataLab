from .utils import load_json, evaluate_typeid, AC_THRESHOLD

GT_PATH = "benchmark/gold_results/themarker_extract_1.json"

def themarker_extract_1(d):
    try:
        PRED_PATH = f"{d}/themarker_extract_1.json"
        gt_lst = load_json(GT_PATH)
        pred_list = load_json(PRED_PATH)
        score = evaluate_typeid(gt_lst, pred_list)
        return score >= AC_THRESHOLD
    except:
        return False