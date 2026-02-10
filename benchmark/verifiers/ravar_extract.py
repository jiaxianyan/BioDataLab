from .utils import load_json, evaluate_ac, AC_THRESHOLD

GT_PATH = "benchmark/gold_results/ravar_extract.json"

def ravar_extract(d):
    try:
        PRED_PATH = f"{d}/ravar_extract.json"
        gt_lst = load_json(GT_PATH)
        pred_list = load_json(PRED_PATH)
        score = evaluate_ac(gt_lst, pred_list, "EFO Trait Label")
        return score >= AC_THRESHOLD
    except:
        return False