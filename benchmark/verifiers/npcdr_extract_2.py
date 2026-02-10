from .utils import evaluate_csv, CVS_THRESHOLD
import pandas as pd

GT_PATH = "benchmark/gold_results/npcdr_extract_2.csv"

def npcdr_extract_2(d):
    try:
        PRED_PATH = f"{d}/npcdr_extract_2.csv"
        gt_df = pd.read_csv(GT_PATH)
        pred_df = pd.read_csv(PRED_PATH)
        score = evaluate_csv(gt_df, pred_df)
        return score >= CVS_THRESHOLD
    except:
        return False