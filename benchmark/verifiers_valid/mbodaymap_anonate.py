import os

def mbodaymap_anonate(d) -> bool:
    results = evaluate_mbodaymap_anonate(
        clean_fasta_path = '../dataset/mBodyMap/SRR1521254_2_clean_R.fasta', 
        raw_fastq_path = '../dataset/mBodyMap/SRR1521254_2.fastq', 
    )
    match = results
    print(f"File exists: {match}")
    return match


def evaluate_mbodaymap_anonate(clean_fasta_path, raw_fastq_path):
    """
    检查生成的文件是否存在。
    
    Args:
        clean_fasta_path (str): 清洗后的 FASTA 文件路径。
        raw_fastq_path (str): 原始 FASTQ 文件路径。
        
    Returns:
        bool: 文件存在返回 True，不存在返回 False。
    """
    return os.path.exists(clean_fasta_path)