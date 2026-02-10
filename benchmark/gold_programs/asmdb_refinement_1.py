import pandas as pd
import os

def filter_asm_candidates(input_path, output_path):
    # 检查文件是否存在
    if not os.path.exists(input_path):
        print(f"错误: 找不到文件 {input_path}")
        return

    try:
        # 读取BED文件，假设没有表头(header=None)，以tab分隔
        # 根据你提供的示例，第11列（索引10）是甲基化百分比（0-100）
        df = pd.read_csv(input_path, sep='\t', header=None)
        
        # 核心逻辑：
        # 1. 获取第11列数据 (df[10])
        # 2. 筛选条件：值 >= 10 (对应0.1) 且 值 <= 90 (对应0.9)
        # 注意：你的示例数据里全是100，这会导致结果为空，这是符合逻辑的
        filtered_df = df[(df[10] >= 10) & (df[10] <= 90)]
        
        # 保存结果
        filtered_df.to_csv(output_path, sep='\t', header=False, index=False)
        
        print(f"处理完成。")
        print(f"原始行数: {len(df)}")
        print(f"保留行数: {len(filtered_df)} (有效ASM候选位点)")
        print(f"结果已保存至: {output_path}")

    except Exception as e:
        print(f"处理过程中发生错误: {e}")

# 定义路径
input_file = "benchmark/dataset/ASMdb/toy_bed/mini_test.bed"
output_file = "benchmark/dataset/ASMdb/toy_bed/mini_test_filtered.bed"

# 执行
filter_asm_candidates(input_file, output_file)