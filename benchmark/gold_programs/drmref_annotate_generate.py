import scanpy as sc
import pandas as pd
import os

def run_benchmark_task():
    # 设置工作目录
    # 1. 加载已经固定好分组的数据集
    # 假设数据就在当前目录下，或者由 benchmark 框架挂载进来
    input_file = "benchmark/dataset/DRMef/pbmc3k_cd4_fixed.h5ad"
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Please make sure {input_file} exists.")
    
    adata = sc.read_h5ad(input_file)
    print(f"Loaded data with {adata.n_obs} cells.")

    # 2. 执行差异表达分析
    # 题目要求：compare 'resistant' against 'sensitive'
    # method: Wilcoxon rank-sum test
    sc.tl.rank_genes_groups(
        adata, 
        groupby='group', 
        method='wilcoxon', 
        groups=['resistant'],   # 感兴趣的组
        reference='sensitive'   # 对照组
    )

    # 3. 提取结果并过滤
    # 获取差异分析的数据框
    result_df = sc.get.rank_genes_groups_df(adata, group='resistant')

    # 过滤条件：
    # adjusted p-value < 0.05
    # absolute log2 fold change > 0.25
    filtered_df = result_df[
        (result_df['pvals_adj'] < 0.05) & 
        (abs(result_df['logfoldchanges']) > 0.25)
    ].copy()

    # 4. 保存结果
    # 要求：single-column CSV, header 'gene'
    output_path = os.path.join('benchmark/dataset/DRMef', "degs.csv")
    
    # scanpy 的结果 DataFrame 中，基因名通常在 'names' 列
    filtered_df[['names']].rename(columns={'names': 'gene'}).to_csv(output_path, index=False)
    
    print(f"Found {len(filtered_df)} genes matching criteria.")
    print(f"Results saved to {output_path}")

if __name__ == "__main__":
    run_benchmark_task()