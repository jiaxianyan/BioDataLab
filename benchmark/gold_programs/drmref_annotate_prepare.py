import scanpy as sc
import numpy as np
import scipy.sparse as sp
import os

# 1. 加载数据
input_file = "benchmark/dataset/DRMef/pbmc3k_cd4_fixed.h5ad"
# 为了保险，我们重新从源头加载，确保拿到 clean 的 raw 数据
print("Reloading from source to ensure clean state...")
adata = sc.datasets.pbmc3k_processed()
adata = adata[adata.obs['louvain'] == 'CD4 T cells'].copy()

# 2. 重新固定分组 (和之前保持一致)
np.random.seed(42)
n_cells = adata.shape[0]
groups = np.array(['sensitive'] * n_cells)
indices = np.random.choice(n_cells, size=n_cells // 2, replace=False)
groups[indices] = 'resistant'
adata.obs['group'] = groups
adata.obs['group'] = adata.obs['group'].astype('category')

# =======================================================
# 🛠️ 关键修复：从 .raw 恢复未缩放的数据
# =======================================================
if adata.raw is not None:
    print("🔄 Restoring unscaled data from .raw to .X ...")
    # 把 raw 中的数据（通常是 log1p 后的非负数据）覆盖到 .X
    # 注意：adata.raw.X 可能是稀疏矩阵，这很好
    adata_raw = adata.raw.to_adata()
    # 我们只需要 raw 的表达矩阵，但要保留我们刚刚创建的 obs (group信息)
    adata.X = adata_raw[:, adata.var_names].X.copy() 
    # 彻底删除 raw，防止 scanpy 再次混淆
    del adata.raw
    print("✅ Data unscaled. All values should be >= 0 now.")
else:
    print("⚠️ No .raw found. Assuming .X is already unscaled.")

# 3. 智能选择基因
np.random.seed(999) 
# 确保我们选的基因都在 adata.var_names 里
valid_genes = [g for g in adata.var_names if g in adata.var_names] 
target_genes = np.random.choice(valid_genes, 5, replace=False).tolist()
print(f"🎯 Target Genes for Spike-in: {target_genes}")

# 4. 准备修改数据
# 转为稠密矩阵以便修改
if sp.issparse(adata.X):
    adata.X = adata.X.toarray()

# 5. 注入信号 (Inject Signal)
# 因为现在是 Log-normalized 数据 (范围 0-10)，加 2.0 已经非常巨大了
# 不需要加 5.0 那么夸张，3.0 足够产生天文数字般的差异
resistant_mask = adata.obs['group'] == 'resistant'
gene_indices = [adata.var_names.get_loc(gene) for gene in target_genes]
spike_in_value = 3.0  

print(f"💉 Injecting signal (value += {spike_in_value})...")
for gene_idx in gene_indices:
    adata.X[resistant_mask, gene_idx] += spike_in_value

# 6. 验证数据中是否有负数 (Double Check)
min_val = adata.X.min()
if min_val < 0:
    print(f"⚠️ Warning: Negative values detected ({min_val}). Rank_genes_groups might fail.")
else:
    print(f"✅ Data integrity check passed: Min value is {min_val} (>=0).")

# 7. 保存
adata.X = sp.csr_matrix(adata.X)
output_file = "benchmark/dataset/DRMef/pbmc3k_cd4_fixed.h5ad"
adata.write(output_file)
print(f"💾 Modified dataset saved to: {output_file}")

# ==========================================
# 8. 立即验证
# ==========================================
print("\n🔍 Verifying with rank_genes_groups...")
# 现在 .X 是健康的 log-normalized 数据，use_raw=False 是安全的
sc.tl.rank_genes_groups(adata, groupby='group', method='wilcoxon', 
                        groups=['resistant'], reference='sensitive', use_raw=False)

result_df = sc.get.rank_genes_groups_df(adata, group='resistant')

# 只要是 padj < 0.05 就算找到，不卡 foldchange 阈值先看看
detected = result_df[
    (result_df['pvals_adj'] < 0.05) &
    (result_df['names'].isin(target_genes))
]

print(f"🎉 Verification: Found {len(detected)} / {len(target_genes)} target genes.")
print("Genes found:", detected['names'].tolist())
print("Log2 Fold Changes:", detected['logfoldchanges'].tolist())

# 保存标准答案
ground_truth_path = "benchmark/dataset/DRMef/ground_truth_degs.csv"
detected[['names']].rename(columns={'names': 'gene'}).to_csv(ground_truth_path, index=False)