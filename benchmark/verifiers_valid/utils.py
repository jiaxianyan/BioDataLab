import os

def evaluate_lst(gt_list, pred_list):
    """检查生成的文件是否存在"""
    return os.path.exists(pred_list)

def evaluate_typeid(gt_list, pred_list, key='PMID'):
    """检查生成的文件是否存在"""
    return os.path.exists(pred_list)

def evaluate_ac(gt_list, pred_list, key='PMID'):
    """检查生成的文件是否存在"""
    return os.path.exists(pred_list)

def evaluate_csv(df_gt, df_pred):
    """检查生成的文件是否存在"""
    return os.path.exists(df_pred)