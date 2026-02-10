from .utils import evaluate_csv, CVS_THRESHOLD
import pandas as pd

GT_PATH = "benchmark/gold_results/dda_extract.csv"

def dda_extract(d):
    try:
        PRED_PATH = f"{d}/dda_extract.csv"
        gt_df = pd.read_csv(GT_PATH)
        pred_df = pd.read_csv(PRED_PATH)
        score = evaluate_csv(gt_df, pred_df)
        return score >= CVS_THRESHOLD
    except:
        return False