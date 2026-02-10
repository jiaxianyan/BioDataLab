import json
import pandas as pd
import numpy as np
from scipy.optimize import linear_sum_assignment

CVS_THRESHOLD = 0.8
AC_THRESHOLD = 1

def load_json(file_path):
    with open(file_path, 'r', encoding='utf-8') as f:
        return json.load(f)

def calculate_metrics(gt_set, pred_set):
    tp = len(gt_set.intersection(pred_set))  # 交集：找对了
    fp = len(pred_set - gt_set)              # 差集：模型多找了/找错了
    fn = len(gt_set - pred_set)              # 差集：模型漏找了

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0.0

    return round(f1, 4)

def evaluate_lst(gt_list, pred_list):
    gt_set = set(gt_list)
    pred_set = set(pred_list)
    
    return calculate_metrics(gt_set, pred_set)

def evaluate_typeid(gt_list, pred_list, key='PMID'):
    def normalize_type(val):
        return str(val).strip()

    gt_map = {
        str(item[key]).strip(): normalize_type(item["TypeID"]) 
        for item in gt_list
    }
    
    pred_map = {
        str(item[key]).strip(): normalize_type(item["TypeID"]) 
        for item in pred_list
    }
    
    total_samples = len(gt_map)
    
    if total_samples == 0:
        return 1.0 if len(pred_map) == 0 else 0.0

    correct_count = 0
    for pmid, true_type in gt_map.items():
        if pmid in pred_map:
            pred_type = pred_map[pmid]
            if pred_type == true_type:
                correct_count += 1
            else:
                pass
        else:
            pass

    accuracy = correct_count / total_samples
    
    return round(accuracy, 4)

def evaluate_ac(gt_list, pred_list, key='PMID'):
    gt_set = {
        item[key] for item in gt_list 
        if item.get('Accepted') is True
    }
    
    pred_set = {
        item[key] for item in pred_list 
        if item.get('Accepted') is True
    }
    
    return calculate_metrics(gt_set, pred_set)

def evaluate_csv(df_gt: pd.DataFrame, df_pred: pd.DataFrame) -> float:
    """
    计算基于匈牙利算法的 Soft F1 Score。
    
    参数:
    df_gt   : Ground Truth DataFrame
    df_pred : Prediction DataFrame
    
    返回:
    float : F1 Score (0.0 ~ 1.0)
    """
    
    df_gt = df_gt.copy()
    df_pred = df_pred.copy()
    df_gt.columns = df_gt.columns.astype(str).str.strip()
    df_pred.columns = df_pred.columns.astype(str).str.strip()

    common_cols = list(set(df_gt.columns).intersection(set(df_pred.columns)))
    
    if not common_cols:
        return 0.0

    # 3. 转换为字典列表 (List of Dicts)，只保留公共列
    #    使用 'records' 模式保留每一行的数据结构
    gt_list = df_gt[common_cols].to_dict('records')
    pred_list = df_pred[common_cols].to_dict('records')

    n_gt = len(gt_list)
    n_pred = len(pred_list)

    if n_gt == 0 and n_pred == 0: return 1.0
    if n_gt == 0 or n_pred == 0: return 0.0

    # ================= 内部辅助函数 =================
    
    def is_cell_match(v1, v2):
        s1 = str(v1).strip()
        s2 = str(v2).strip()
        
        if s1.startswith('"') and s1.endswith('"'):
            s1_unquoted = s1[1:-1]
            return 1.0 if s1_unquoted.lower() == s2.lower() else 0.0

        if ';' in s1:
            synonyms = [syn.strip().lower() for syn in s1.split(';')]
            for syn in synonyms:
                if s2.lower() == syn:
                    return 1.0
            return 0.0
        else:
            return 1.0 if s1.lower() == s2.lower() else 0.0

    def calculate_row_similarity(row_a, row_b):
        matches = 0.0
        for col in common_cols:
            matches += is_cell_match(row_a.get(col), row_b.get(col))
        return matches / len(common_cols)

    # ================= 核心计算逻辑 =================
    
    # 5. 构建相似度矩阵 [GT行数 x Pred行数]
    sim_matrix = np.zeros((n_gt, n_pred))

    for i in range(n_gt):
        for j in range(n_pred):
            sim_matrix[i, j] = calculate_row_similarity(gt_list[i], pred_list[j])

    # 6. 匈牙利算法求解最大权匹配
    # linear_sum_assignment 寻找最小成本，因此 Cost = 1.0 - Similarity
    cost_matrix = 1.0 - sim_matrix
    row_ind, col_ind = linear_sum_assignment(cost_matrix)

    # 7. 计算指标
    # 提取最佳匹配的总相似度得分
    total_match_score = sim_matrix[row_ind, col_ind].sum()

    # Precision: (匹配总分) / (预测行数) -> 惩罚多预测的行
    precision = total_match_score / n_pred
    
    # Recall: (匹配总分) / (GT行数) -> 惩罚漏预测的行
    recall = total_match_score / n_gt

    if (precision + recall) == 0:
        return 0.0

    f1 = 2 * (precision * recall) / (precision + recall)
    
    return round(f1, 4)


