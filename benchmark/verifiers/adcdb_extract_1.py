from .utils import evaluate_csv, CVS_THRESHOLD
import pandas as pd

GT_PATH = "benchmark/gold_results/adcdb_extract_1.csv"


def adcdb_extract_1(d):
    try:
        PRED_PATH = f"{d}/adcdb_extract_1.csv"
        gt_df = pd.read_csv(GT_PATH)
        pred_df = pd.read_csv(PRED_PATH)
        score = evaluate_csv(gt_df, pred_df)
        return score >= CVS_THRESHOLD
    except:
        return False