import json

def cancerscem_annotate(d):
    try:
        # 读取模型预测结果
        with open(f'{d}/cancerscem_annotate.json', 'r') as f:
            pred_result = json.load(f)

        # 读取标准答案
        with open('benchmark/gold_results/cancerscem_annotate.json', 'r') as f:
            gold_result = json.load(f)

        # 直接比较两个 dict，严格区分大小写
        return pred_result == gold_result

    except Exception as e:
        print(f"Error processing cancerscem annotation results: {e}")
        return False
