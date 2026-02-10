import os
import pandas as pd
import re

def metazexp_annotate(d):
    """
    Evaluates DESeq2 results by comparing the overlap of top significant genes (Method 2).
    
    Args:
        d (str): Directory containing the predicted 'deseq2_results.tsv'.
        
    Returns:
        bool: True if the top significant genes match the ground truth with high overlap.
    """
    # ================= 配置区域 =================
    PRED_FILENAME = 'metazexp_annotate.tsv'       # 预测结果文件名
    GOLD_PATH = 'benchmark/gold_results/metazexp_annotate.tsv' # 金标准路径
    TOP_N = 50                                 # 只比较最显著的前 50 个基因
    JACCARD_THRESHOLD = 0.9                    # 重合度阈值 (0.9 表示 90% 重合即算对)
    # ===========================================

    pred_path = os.path.join(d, PRED_FILENAME)

    try:
        # 1. 检查文件是否存在
        if not os.path.exists(pred_path):
            print(f"[Fail] Prediction file not found: {pred_path}")
            return False

        # 2. 读取数据 (假设是 Tab 分隔)
        df_pred = pd.read_csv(pred_path, sep='\t')
        df_gold = pd.read_csv(GOLD_PATH, sep='\t')

        # 3. 关键列检查
        required_cols = ['gene_id', 'padj']
        if not all(col in df_pred.columns for col in required_cols):
            print(f"[Fail] Missing columns. Required: {required_cols}")
            return False

        # 4. 数据清洗 (处理 'e-64x' 这种脏数据)
        def clean_numeric(series):
            # 将列转换为字符串，移除所有非数字、非小数点、非科学计数法符号(e/E/-)的字符
            # 然后强制转换为 float，无法转换的变为 NaN
            clean_str = series.astype(str).astype(str).str.replace(r'[^\d\.eE-]', '', regex=True)
            return pd.to_numeric(clean_str, errors='coerce')

        # 对 padj 列进行清洗
        df_pred['padj'] = clean_numeric(df_pred['padj'])
        df_gold['padj'] = clean_numeric(df_gold['padj'])

        # 5. 提取 Top-N 显著基因
        # 逻辑：去除 NaN 值 -> 按 padj 从小到大排序 -> 取前 N 个 -> 获取 gene_id
        def get_top_genes(df, n):
            # dropna: 移除没有计算出 p-value 的行
            sorted_df = df.dropna(subset=['padj']).sort_values(by='padj', ascending=True)
            # 取前 n 个基因 ID，转换为集合
            return set(sorted_df.head(n)['gene_id'].astype(str).str.strip().str.lower())

        gold_genes_set = get_top_genes(df_gold, TOP_N)
        pred_genes_set = get_top_genes(df_pred, TOP_N)

        # 6. 计算 Jaccard Index (交集 / 并集)
        # 或者使用简单的 Recall (交集 / 金标准数量)
        intersection = len(gold_genes_set.intersection(pred_genes_set))
        union = len(gold_genes_set.union(pred_genes_set))
        
        if union == 0:
            print("[Fail] No valid genes found in both lists.")
            return False

        jaccard_index = intersection / union
        
        # 打印详细信息方便调试
        print(f"Top-{TOP_N} Evaluation:")
        print(f"  - Intersection: {intersection}")
        print(f"  - Union: {union}")
        print(f"  - Jaccard Index: {jaccard_index:.4f} (Threshold: {JACCARD_THRESHOLD})")

        # 7. 判定结果
        if jaccard_index >= JACCARD_THRESHOLD:
            return True
        else:
            return False

    except Exception as e:
        print(f"[Error] An error occurred during evaluation: {e}")
        return False