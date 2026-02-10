import pandas as pd

# 1. 读取数据
print("1. 读取数据...")
genes = pd.read_csv("benchmark/dataset/ASMdb/mock/mock_genes.csv")
asm = pd.read_csv("benchmark/dataset/ASMdb/mock/mock_asm_sites.csv")

# 2. 计算每个基因的搜索窗口 (Gene Body + Upstream 3kb)
# 逻辑：
# 正链 (+): [Start - 3000, End]
# 负链 (-): [Start, End + 3000]
print("2. 计算搜索窗口...")
genes['search_start'] = genes.apply(
    lambda x: x['start'] - 3000 if x['strand'] == '+' else x['start'], axis=1
)
genes['search_end'] = genes.apply(
    lambda x: x['end'] if x['strand'] == '+' else x['end'] + 3000, axis=1
)
# 确保坐标不小于0
genes['search_start'] = genes['search_start'].clip(lower=0)

# 3. 空间重叠判定 (Spatial Join)
# 这是一个简化的 Python 循环实现。大数据量通常使用 bedtools 或 pybedtools
print("3. 进行空间映射 (这可能需要几秒钟)...")

# 将 ASM 数据按染色体分组，加速查询
asm_grouped = asm.groupby('chrom')

results = []

for idx, gene in genes.iterrows():
    chrom = gene['chrom']
    
    # 如果该染色体上有 ASM 位点
    if chrom in asm_grouped.groups:
        chrom_asm = asm_grouped.get_group(chrom)
        
        # 核心逻辑：找出落在 [search_start, search_end] 范围内的位点
        # 使用 Series 的 between 方法
        hits = chrom_asm[chrom_asm['pos'].between(gene['search_start'], gene['search_end'])]
        
        if not hits.empty:
            # 统计有多少个 *不同* 的样本在这个基因附近有 ASM
            # 原文需求："将基因覆盖ASM的频率在所有样本中相加" -> 即 Count Unique Samples
            unique_samples = hits['sample_id'].nunique()
            results.append({'gene_id': gene['gene_id'], 'frequency': unique_samples})
        else:
            results.append({'gene_id': gene['gene_id'], 'frequency': 0})
    else:
        results.append({'gene_id': gene['gene_id'], 'frequency': 0})

df_results = pd.DataFrame(results)

# 4. 排序并取 Top 100
print("4. 排序并输出结果...")
top_100 = df_results.sort_values(by='frequency', ascending=False).head(100)

print("\n=== Top 10 ASM Genes ===")
print(top_100.head(10))

# 导出结果
top_100.to_csv("benchmark/dataset/ASMdb/mock/benchmark_result_top100.csv", index=False)
print("\n✅ 结果已保存至 benchmark_result_top100.csv")