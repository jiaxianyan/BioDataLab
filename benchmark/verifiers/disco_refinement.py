import os
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr

def disco_refinement(d):
    """
    Evaluates whether the generated RDS file matches the ground truth.
    
    Args:
        d (str): The directory containing the generated 'disco_refinement.rds'.
                 Example: 'benchmark/output/task_1'
                 
    Returns:
        bool: True if files are semantically identical, False otherwise.
    """
    try:
        # 1. Construct Paths
        # 生成结果路径 (根据题目描述，结果在 <r></r>/disco_refinement.rds，这里 d 代表该路径)
        pred_path = os.path.join(d, 'disco_refinement.rds')
        # 金标文件路径 (假设位于 standard benchmark location)
        gold_path = 'benchmark/gold_results/disco_refinement.rds'

        # 2. Check File Existence
        if not os.path.exists(pred_path):
            print(f"[Evaluation] Generated file not found: {pred_path}")
            return False

        # 3. Load R Interface
        # 使用 R 的 base 包来读取和比较
        base = importr('base')
        
        # 4. Read RDS Files
        # readRDS 可以正确加载 Matrix 包定义的 S4 稀疏矩阵
        try:
            pred_obj = base.readRDS(pred_path)
            gold_obj = base.readRDS(gold_path)
        except Exception as e:
            print(f"[Evaluation] Error reading RDS files (R format error): {e}")
            return False

        # 5. Compare Objects
        # 使用 R 原生的 all.equal() 函数。
        # 它非常强大，能够自动处理：
        #   - 稀疏矩阵的数值比较
        #   - 行列名 (Dimnames) 的一致性
        #   - 维度的匹配
        #   - 能够容忍浮点数的微小误差 (默认 tolerance)
        # 注意：如果格式一个是 dgTMatrix 一个是 dgCMatrix 但内容一致，all.equal 通常会视为相等(取决于R版本和包)，
        # 如需严格限制格式，可增加 class 检查。
        
        comparison_result = base.all_equal(pred_obj, gold_obj)

        # all.equal 返回 True (bool) 或者 描述差异的字符串列表
        # 在 rpy2 中，True 会被转换为包含 True 的 R 向量
        
        # 将 R 向量转换为 Python 元组/列表进行判断
        res_py = tuple(comparison_result)
        
        # 只有当结果仅包含 Python 的 True 时，才算通过
        if len(res_py) == 1 and res_py[0] is True:
            return True
        else:
            # 打印差异信息以便调试 (可选)
            # print(f"[Evaluation] Mismatch details: {res_py}") 
            return False

    except Exception as e:
        print(f"[Evaluation] System error: {e}")
        return False