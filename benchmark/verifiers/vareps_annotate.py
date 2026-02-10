
def vareps_annotation(d):
    try:
        # 预测结果路径
        pred_path = f"{d}/most_stable_mutation.txt"
        # 标准答案路径
        gold_path = "benchmark/gold_results/vareps_annotate.txt"

        with open(pred_path, "r") as f:
            pred_text = f.read().strip()

        with open(gold_path, "r") as f:
            gold_text = f.read().strip()

        # 只判断文本是否一致
        return pred_text == gold_text

    except Exception as e:
        print(f"Error in vareps_annotation: {e}")
        return False
