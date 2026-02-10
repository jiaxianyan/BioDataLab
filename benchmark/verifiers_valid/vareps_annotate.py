import os

def vareps_annotation(d):
    try:
        # 预测结果路径
        pred_path = f"{d}/most_stable_mutation.txt"
        
        # 检查文件是否存在
        return os.path.exists(pred_path)
    
    except Exception as e:
        print(f"Error in vareps_annotation: {e}")
        return False