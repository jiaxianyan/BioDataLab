from .utils import load_json, evaluate_ac, AC_THRESHOLD

GT_PATH = "benchmark/gold_results/adcdb_extract_2.json"
# PRED_PATH = "<r></r>/adcdb_extract_2.json"

def adcdb_extract_2(d):
    try:
        PRED_PATH = f"{d}/adcdb_extract_2.json"
        gt_lst = load_json(GT_PATH)
        pred_list = load_json(PRED_PATH)
        score = evaluate_ac(gt_lst, pred_list)
        return score >= AC_THRESHOLD
    except:
        return False