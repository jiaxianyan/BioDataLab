import gzip
import os

def mbodaymap_anonate(d) -> bool:
    results = evaluate_mbodaymap_anonate(
        clean_fasta_path = '../dataset/mBodyMap/SRR1521254_2_clean_R.fasta', 
        raw_fastq_path = '../dataset/mBodyMap/SRR1521254_2.fastq', 
    )
    if results["is_compliant"]==1:
         match = True
    else:
        match = False
    print(results)
    return match


def evaluate_mbodaymap_anonate(clean_fasta_path, raw_fastq_path):
    """
    评估清洗后的 FASTA 数据是否符合文献要求的质量标准。
    
    Args:
        clean_fasta_path (str): 清洗后的 FASTA 文件路径。
        raw_fastq_path (str): 原始 FASTQ 文件路径 (支持 .gz)，用于确定原始读长。
        
    Returns:
        dict: 包含数值化评估指标的字典 metrics。
    """
    metrics = {
        "total_reads": 0,       # 总序列数
        "threshold_bp": 0,      # 计算出的 2/3 长度阈值
        "is_empty": 1,          # 1=空文件, 0=非空 (数值指标)
        "is_compliant": 0,      # 1=全部通过, 0=存在失败 (数值指标)
        "failure_count": 0,     # 低于阈值的序列数
        "failure_rate": 0.0     # 低于阈值的序列占比 (0.0 - 1.0)
    }

    # 1. 获取原始读长并计算阈值
    # -------------------------------------------------
    if not os.path.exists(raw_fastq_path):
        raise FileNotFoundError(f"Raw file not found: {raw_fastq_path}")
    
    try:
        # 兼容普通文件和 gzip 文件
        opener = gzip.open if raw_fastq_path.endswith('.gz') else open
        with opener(raw_fastq_path, 'rt') as f:
            # FASTQ 格式：第2行是序列
            _ = f.readline() # Skip ID
            seq_line = f.readline()
            if not seq_line:
                raise ValueError("Raw FASTQ file is empty or invalid.")
            original_len = len(seq_line.strip())
            
            # 文献要求: "Sequences shorter than two-thirds... excluded"
            # 使用 int() 向下取整
            threshold = int(original_length := original_len * 2 / 3)
            metrics["threshold_bp"] = threshold
    except Exception as e:
        print(f"Error reading raw file: {e}")
        return metrics

    # 2. 扫描 Clean FASTA 文件并统计
    # -------------------------------------------------
    if not os.path.exists(clean_fasta_path):
        # 如果文件不存在，保持 is_empty=1
        return metrics

    total_seqs = 0
    failed_seqs = 0
    
    # 简单的状态机解析 FASTA (处理多行序列的情况)
    current_seq_len = 0
    in_sequence = False

    try:
        # 同样兼容 .gz (虽然一般 clean data 可能是解压的，这里做个防御)
        opener = gzip.open if clean_fasta_path.endswith('.gz') else open
        with opener(clean_fasta_path, 'rt') as f:
            for line in f:
                line = line.strip()
                if not line: continue
                
                if line.startswith(">"):
                    # 结算上一条序列
                    if in_sequence:
                        total_seqs += 1
                        if current_seq_len < threshold:
                            failed_seqs += 1
                    
                    # 重置下一条
                    current_seq_len = 0
                    in_sequence = True
                else:
                    # 累加序列长度
                    current_seq_len += len(line)
            
            # 结算最后一条序列
            if in_sequence:
                total_seqs += 1
                if current_seq_len < threshold:
                    failed_seqs += 1
                    
    except Exception as e:
        print(f"Error reading fasta file: {e}")
        return metrics
    
    metrics["total_reads"] = total_seqs
    metrics["failure_count"] = failed_seqs
    
    if total_seqs > 0:
        metrics["is_empty"] = 0 # 非空
        metrics["failure_rate"] = round(failed_seqs / total_seqs, 6)
        
        # 判断是否合规：要求 failed_seqs 必须为 0
        if failed_seqs == 0:
            metrics["is_compliant"] = 1
        else:
            metrics["is_compliant"] = 0
    else:
        metrics["is_empty"] = 1 # 空文件
        metrics["is_compliant"] = 0 # 空文件视为不合规（或根据需求调整）

    return metrics