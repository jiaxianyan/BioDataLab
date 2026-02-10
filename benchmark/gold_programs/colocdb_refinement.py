import pandas as pd
import os
import sys

# ================= 配置部分 =================
input_file = "GCST90000064_buildGRCh37.tsv"
output_file = "GCST90000064_buildGRCh38.tsv"
chain_file = "hg19ToHg38.over.chain.gz"
liftover_binary = "liftOver"  # 确保 liftover 有执行权限 (chmod +x liftover)

# 临时文件命名
temp_bed_input = "temp_input.bed"
temp_bed_output = "temp_output.bed"
temp_bed_unmapped = "temp_unmapped.bed"

def main():
    print(f"1. 正在读取原始文件: {input_file} ...")
    df = pd.read_csv(input_file, sep='\t')
    
    # 检查必要的列是否存在
    required_cols = ['chromosome', 'base_pair_location']
    if not all(col in df.columns for col in required_cols):
        print(f"错误: 输入文件缺少必要的列: {required_cols}")
        return

    print("2. 正在生成临时 BED 文件用于 liftOver ...")
    # 创建用于 liftover 的 DataFrame
    # 格式: chr, start, end, id (我们使用原始 dataframe 的 index 作为 id 以便后续合并)
    bed_df = pd.DataFrame()
    
    # 处理染色体: 添加 'chr' 前缀，并转为字符串 (例如 1 -> chr1)
    # 注意: 如果你的 chain 文件不需要 chr 前缀，请去除 'chr' + 
    bed_df['chrom'] = 'chr' + df['chromosome'].astype(str)
    
    # 处理坐标: BED 是 0-based, GWAS 是 1-based
    # 转换: start = pos - 1, end = pos
    bed_df['start'] = df['base_pair_location'] - 1
    bed_df['end'] = df['base_pair_location']
    
    # 添加 ID (行索引) 用来把转换后的结果对应回去
    bed_df['id'] = df.index
    
    # 保存为无表头的 BED 文件
    bed_df.to_csv(temp_bed_input, sep='\t', header=False, index=False)
    
    print("3. 正在运行 liftOver 工具 ...")
    # 构建命令
    cmd = f"{liftover_binary} {temp_bed_input} {chain_file} {temp_bed_output} {temp_bed_unmapped}"
    
    # 执行命令
    exit_code = os.system(cmd)
    if exit_code != 0:
        print("错误: liftOver 运行失败，请检查工具路径和 chain 文件是否正确。")
        return

    print("4. 正在处理转换结果并合并数据 ...")
    # 读取转换成功的 BED 文件
    # 输出格式: chrom, start, end, id
    try:
        mapped_bed = pd.read_csv(temp_bed_output, sep='\t', header=None, names=['chrom', 'start', 'end', 'id'])
    except pd.errors.EmptyDataError:
        print("错误: 没有数据成功转换 (输出文件为空)。请检查 Chain 文件是否匹配。")
        return

    # 统计转换率
    total_count = len(df)
    mapped_count = len(mapped_bed)
    print(f"   转换统计: {mapped_count}/{total_count} 个位点成功转换 ({mapped_count/total_count:.2%})")

    # 将转换后的坐标映射回原始数据
    # 我们使用 'id' (原索引) 进行合并
    # 转换后的 GWAS 坐标应该等于 BED 的 'end' 列 (变回 1-based)
    
    # 创建一个字典: id -> new_pos
    id_to_pos = dict(zip(mapped_bed['id'], mapped_bed['end']))
    
    # 筛选出成功的行
    df_clean = df.loc[mapped_bed['id']].copy()
    
    # 更新坐标
    df_clean['base_pair_location'] = df_clean.index.map(id_to_pos)
    
    # 可选：如果你希望输出文件保持不带 'chr' 的格式 (如 '1' 而不是 'chr1')，不需要改动 chromosome 列
    # 因为位置变了，但染色体编号通常对应的还是原来的数字。
    # 除非发生了染色体间的重排，通常我们信任原始的 chromosome 列，或者从 mapped_bed 更新它。
    # 严谨做法是从 mapped_bed 更新 chromosome，并去掉 'chr' 前缀:
    id_to_chrom = dict(zip(mapped_bed['id'], mapped_bed['chrom']))
    new_chroms = df_clean.index.map(id_to_chrom).str.replace('chr', '')
    df_clean['chromosome'] = new_chroms

    print(f"5. 正在保存结果到: {output_file} ...")
    df_clean.to_csv(output_file, sep='\t', index=False)

    # 清理临时文件
    print("6. 清理临时文件 ...")
    # os.remove(temp_bed_input)
    # os.remove(temp_bed_output)
    # temp_unmapped.bed 可以保留供你检查哪些位点丢失了
    # os.remove(temp_bed_unmapped) 
    
    print("完成！")

if __name__ == "__main__":
    main()