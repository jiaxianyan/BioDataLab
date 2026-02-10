import json

def mbodymap_integration(d):
    try:
        # 读取模型生成的结果
        with open(f'{d}/mbodymap_integration.json', 'r') as f:
            pred_srr_list = json.load(f)

        # 读取 benchmark 中的标准答案
        with open('benchmark/gold_results/mbodymap_integration.json', 'r') as f:
            gold_srr_list = json.load(f)

        # 不区分顺序，区分大小写，直接排序比较
        return sorted(pred_srr_list) == sorted(gold_srr_list)

    except Exception as e:
        print(f"Error processing JSON files: {e}")
        return False
