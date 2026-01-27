import json
import hashlib
import re
import sys
import os
from typing import Dict, Optional, Tuple
import subprocess
import math
import pandas as pd
import numpy as np

def eval_str_from_file_equal(file1_path: str, file2_path: str) -> bool:
    """
    Evaluates if two files contain the exact same content using the filecmp module.

    This is the recommended, robust, and Pythonic way to compare files.

    Args:
        file1_path: The path to the first file.
        file2_path: The path to the second file.

    Returns:
        True if the files have identical content, False otherwise.
        Returns False if one or both files do not exist.
    """
    # First, check if both files exist to avoid filecmp raising an error
    if not os.path.exists(file1_path) or not os.path.exists(file2_path):
        print("Error: One or both files do not exist.")
        return 

    with open(file1_path, 'r') as f:
        text1 = f.read().strip()

    with open(file2_path, 'r') as f:
        text2 = f.read().strip()

    # shallow=False ensures a full byte-by-byte content comparison
    return text1 == text2


def get_flagstat_info(filepath: str) -> Optional[Dict[str, Dict[str, int]]]:
    """
    Runs `samtools flagstat` on a given file and parses the output.

    Args:
        filepath (str): The path to the SAM or BAM file.

    Returns:
        A dictionary containing the parsed flagstat information, or None if an
        error occurs. The structure is:
        {
            "in total": {"passed": X, "failed": Y},
            "secondary": {"passed": X, "failed": Y},
            ...
        }
    """
    command = ['samtools', 'flagstat', filepath]
    stats_dict = {}

    try:
        # Execute the samtools command
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=True,  # This will raise CalledProcessError on non-zero exit codes
            encoding='utf-8'
        )
    except FileNotFoundError:
        print(f"Error: 'samtools' command not found.", file=sys.stderr)
        print("Please ensure samtools is installed and in your system's PATH.", file=sys.stderr)
        return None
    except subprocess.CalledProcessError as e:
        print(f"Error running samtools on '{filepath}':", file=sys.stderr)
        print(e.stderr, file=sys.stderr)
        return None

    # Regex to capture "passed + failed" numbers and the description
    # e.g., "3020056 + 0 in total (QC-passed reads + QC-failed reads)"
    pattern = re.compile(r"(\d+) \+ (\d+) (.+)")

    for line in result.stdout.strip().split('\n'):
        match = pattern.match(line)
        if match:
            passed_reads = int(match.group(1))
            failed_reads = int(match.group(2))
            
            # Clean up the description key
            description = match.group(3).strip()
            # Remove extra details in parentheses
            description_key = description.split('(')[0].strip()

            stats_dict[description_key] = {
                "passed": passed_reads,
                "failed": failed_reads
            }
            
    if not stats_dict:
        print(f"Warning: Could not parse any flagstat output for '{filepath}'.", file=sys.stderr)
        return None
        
    return stats_dict

def eval_sam_bam_flagstat(file1_path: str, file2_path: str) -> bool:
    """Main function to parse arguments and compare files."""

    print(f"Comparing '{file1_path}' and '{file2_path}'...\n")

    stats1 = get_flagstat_info(file1_path)
    stats2 = get_flagstat_info(file2_path)

    # Exit if there was an error getting stats for either file
    if not stats1 or not stats2:
        return False

    # Direct comparison of the parsed dictionaries
    if stats1 == stats2:
        print("Success: The samtools flagstat outputs are identical.")
        return True
    else:
        return False

def eval_file_exist(file1_path: str, file2_path: str):
    if os.path.exists(file1_path) and os.path.exists(file2_path):
        return True
    else:
        return False

def eval_numeric_value_from_file_equal_round(file_path1: str, file_path2: str, abs_tol: float = 1e-6) -> bool:
    """
    从两个文件中分别读取一个数值，并判断它们是否在指定的绝对容差内近似相等。

    Args:
        file_path1 (str): 第一个文件的路径。
        file_path2 (str): 第二个文件的路径。
        abs_tol (float, optional): 用于比较的绝对容差。默认为 1e-9。
                                   如果两个数字的差的绝对值小于或等于此值，则认为它们相等。

    Returns:
        bool: 如果文件存在、内容为有效数字且两个数字近似相等，则返回 True。
              在任何错误（如文件未找到、内容无法转换为数字）或数字不相等的情况下，返回 False。
    """
    try:
        # --- 读取第一个文件 ---
        with open(file_path1, 'r') as f1:
            # .strip() 用于移除数字周围可能存在的空格或换行符
            content1 = f1.read().strip()
            num1 = float(content1)

        # --- 读取第二个文件 ---
        with open(file_path2, 'r') as f2:
            content2 = f2.read().strip()
            num2 = float(content2)

        # --- 使用 math.isclose() 进行比较 ---
        # 这是处理浮点数比较的标准方法
        return math.isclose(num1, num2, abs_tol=abs_tol)

    except FileNotFoundError as e:
        print(f"错误：文件未找到 -> {e.filename}")
        return False
    except ValueError:
        # 如果文件内容不是有效的数字（例如包含文本），float()会抛出ValueError
        print(f"错误：文件内容无法转换为数字。请确保文件只包含一个数字。")
        return False
    except Exception as e:
        # 捕获其他可能的异常，例如权限错误
        print(f"发生未知错误: {e}")
        return False

def eval_str_list_from_json_equal_unsort(file1_path: str, file2_path: str) -> bool:
    """
    从两个JSON文件中加载字符串列表，并比较它们是否完全相等（顺序和内容都必须相同）。

    Args:
        file1_path: 第一个JSON文件的路径，应包含一个字符串列表。
        file2_path: 第二个JSON文件的路径，应包含一个字符串列表。

    Returns:
        如果两个列表在顺序和内容上都完全相同，则返回 True，否则返回 False。
        如果文件不存在、不是有效的JSON或内容不是列表，则返回 False。
    """
    if not os.path.exists(file1_path) or not os.path.exists(file2_path):
        print(f"Error: One or both files do not exist ('{file1_path}', '{file2_path}').")
        return False

    try:
        with open(file1_path, 'r', encoding='utf-8') as f:
            list1: List[str] = json.load(f)

        with open(file2_path, 'r', encoding='utf-8') as f:
            list2: List[str] = json.load(f)

        # 确保加载的数据是列表类型
        if not isinstance(list1, list) or not isinstance(list2, list):
            print("Error: One or both JSON files do not contain a list at the root.")
            return False
            
        # 直接比较列表，会检查顺序和内容
        return list1 == list2
    except (json.JSONDecodeError, IOError) as e:
        print(f"Error processing JSON files: {e}")
        return False

def eval_str_list_from_json_equal_sort(file1_path: str, file2_path: str) -> bool:
    """
    从两个JSON文件中加载字符串列表，并比较它们的内容是否相等（忽略顺序）。

    Args:
        file1_path: 第一个JSON文件的路径，应包含一个字符串列表。
        file2_path: 第二个JSON文件的路径，应包含一个字符串列表。

    Returns:
        如果两个列表包含完全相同的元素（无论顺序如何），则返回 True，否则返回 False。
        如果文件不存在、不是有效的JSON或内容不是列表，则返回 False。
    """
    if not os.path.exists(file1_path) or not os.path.exists(file2_path):
        print(f"Error: One or both files do not exist ('{file1_path}', '{file2_path}').")
        return False

    try:
        with open(file1_path, 'r', encoding='utf-8') as f:
            list1: List[str] = json.load(f)

        with open(file2_path, 'r', encoding='utf-8') as f:
            list2: List[str] = json.load(f)
            
        if not isinstance(list1, list) or not isinstance(list2, list):
            print("Error: One or both JSON files do not contain a list at the root.")
            return False

        # 先对列表进行排序，然后比较，这样可以忽略原始顺序
        return sorted(list1) == sorted(list2)
    except Exception as e:
        print(f"Error processing JSON files: {e}")
        return False

def eval_numeric_value_list_from_file_equal_round(file1_path: str, file2_path: str, precision: int = 6) -> bool:
    """
    从两个文本文件中加载数值列表（每行一个数字），并在指定的精度上进行四舍五入后进行比较。

    Args:
        file1_path: 第一个文件的路径，每行应包含一个数字。
        file2_path: 第二个文件的路径，每行应包含一个数字。
        precision: 用于四舍五入的小数位数，默认为6。

    Returns:
        如果两个文件包含相同数量的数字，并且对应位置的数字在四舍五入后相等，则返回 True。
        否则返回 False。如果文件不存在或包含非数值内容，也返回 False。
    """
    if not os.path.exists(file1_path) or not os.path.exists(file2_path):
        print(f"Error: One or both files do not exist ('{file1_path}', '{file2_path}').")
        return False

    try:
        with open(file1_path, 'r', encoding='utf-8') as f:
            # 过滤掉空行
            lines1 = [line.strip() for line in f if line.strip()]
            nums1 = [float(line) for line in lines1]

        with open(file2_path, 'r', encoding='utf-8') as f:
            lines2 = [line.strip() for line in f if line.strip()]
            nums2 = [float(line) for line in lines2]

        # 如果数字数量不同，则列表不相等
        if len(nums1) != len(nums2):
            return False

        # 逐个比较四舍五入后的值
        for n1, n2 in zip(nums1, nums2):
            if round(n1, precision) != round(n2, precision):
                return False
        
        # 所有元素都通过了检查
        return True
    except Exception as e:
        print(f"Error processing numeric files: {e}")
        return False

def eval_dict_from_json_equal(file1_path: str, file2_path: str) -> bool:
    """
    从两个JSON文件中加载字典，并比较它们的内容是否完全相等。

    Python的字典比较本身就是递归的，并且不关心键的顺序。

    Args:
        file1_path: 第一个JSON文件的路径，应包含一个字典。
        file2_path: 第二个JSON文件的路径，应包含一个字典。

    Returns:
        如果两个字典具有相同的键和对应的值，则返回 True，否则返回 False。
        如果文件不存在、不是有效的JSON或内容不是字典，则返回 False。
    """
    if not os.path.exists(file1_path) or not os.path.exists(file2_path):
        print(f"Error: One or both files do not exist ('{file1_path}', '{file2_path}').")
        return False

    try:
        with open(file1_path, 'r', encoding='utf-8') as f:
            dict1: Dict[str, Any] = json.load(f)

        with open(file2_path, 'r', encoding='utf-8') as f:
            dict2: Dict[str, Any] = json.load(f)

        if not isinstance(dict1, dict) or not isinstance(dict2, dict):
            print("Error: One or both JSON files do not contain a dictionary at the root.")
            return False

        #直接比较字典。Python会自动处理键的顺序和值的递归比较。
        return dict1 == dict2
    except Exception as e:
        print(f"Error processing JSON files: {e}")
        return False


def eval_csv_file_equal_sort(file1_path, file2_path, ignore_header=False, ignore_row_order=False, ignore_column_order=True):
    """
    使用 pandas 库比较两个CSV文件的内容是否相等。功能更强大。

    Args:
        file1_path: 第一个CSV文件的路径。
        file2_path: 第二个CSV文件的路径。
        ignore_header: 如果为 True，则将第一行视为数据，而不是列名。
        ignore_row_order: 如果为 True，则忽略数据行的顺序进行比较。
        ignore_column_order: 如果为 True，则忽略列的顺序进行比较。
        **kwargs: 传递给 pandas.read_csv 的其他参数，例如 sep=';'。

    Returns:
        如果文件内容根据指定规则相等，则返回 True，否则返回 False。
        如果文件不存在、pandas未安装或解析错误，则返回 False。
    """
    if pd is None:
        return False # pandas 未安装
        
    if not os.path.exists(file1_path) or not os.path.exists(file2_path):
        print(f"Error: One or both files do not exist ('{file1_path}', '{file2_path}').")
        return False

    try:
        # 如果忽略表头，read_csv 应使用 header=None
        header_option = None if ignore_header else 'infer'
        
        df1 = pd.read_csv(file1_path, header=header_option)
        df2 = pd.read_csv(file2_path, header=header_option)

        # 忽略列序
        if ignore_column_order:
            # 检查列名集合是否相同
            if set(df1.columns) != set(df2.columns):
                return False
            # 将 df2 的列序调整为与 df1 一致
            df2 = df2[df1.columns]

        # 忽略行序
        if ignore_row_order:
            # 根据所有列对 DataFrame 进行排序，并重置索引
            df1 = df1.sort_values(by=list(df1.columns)).reset_index(drop=True)
            df2 = df2.sort_values(by=list(df2.columns)).reset_index(drop=True)
            
        # 使用 pandas 内置的 equals 方法进行精确比较
        return df1.equals(df2)
        
    except Exception as e:
        print(f"An error occurred with pandas: {e}")
        return False

def eval_csv_file_equal_unsort(file1_path, file2_path, ignore_header=False, ignore_row_order=True, ignore_column_order=True):
    """
    使用 pandas 库比较两个CSV文件的内容是否相等。功能更强大。

    Args:
        file1_path: 第一个CSV文件的路径。
        file2_path: 第二个CSV文件的路径。
        ignore_header: 如果为 True，则将第一行视为数据，而不是列名。
        ignore_row_order: 如果为 True，则忽略数据行的顺序进行比较。
        ignore_column_order: 如果为 True，则忽略列的顺序进行比较。

    Returns:
        如果文件内容根据指定规则相等，则返回 True，否则返回 False。
        如果文件不存在、pandas未安装或解析错误，则返回 False。
    """
    if pd is None:
        return False # pandas 未安装
        
    if not os.path.exists(file1_path) or not os.path.exists(file2_path):
        print(f"Error: One or both files do not exist ('{file1_path}', '{file2_path}').")
        return False

    try:
        # 如果忽略表头，read_csv 应使用 header=None
        header_option = None if ignore_header else 'infer'
        
        df1 = pd.read_csv(file1_path, header=header_option)
        df2 = pd.read_csv(file2_path, header=header_option)

        # 忽略列序
        if ignore_column_order:
            # 检查列名集合是否相同
            if set(df1.columns) != set(df2.columns):
                return False
            # 将 df2 的列序调整为与 df1 一致
            df2 = df2[df1.columns]

        # 忽略行序
        if ignore_row_order:
            # 根据所有列对 DataFrame 进行排序，并重置索引
            df1 = df1.sort_values(by=list(df1.columns)).reset_index(drop=True)
            df2 = df2.sort_values(by=list(df2.columns)).reset_index(drop=True)
            
        # 使用 pandas 内置的 equals 方法进行精确比较
        return df1.equals(df2)
        
    except Exception as e:
        print(f"An error occurred with pandas: {e}")
        return False

def parse_fasta(file_path):
    """
    Parses a FASTA file and returns a dictionary.

    Args:
        file_path (str): The path to the FASTA file.

    Returns:
        dict: A dictionary where keys are sequence headers (without '>')
              and values are the corresponding sequences.
    """
    sequences = {}
    try:
        with open(file_path, 'r') as f:
            header = None
            current_sequence = []
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # If we have a previous sequence, save it
                    if header:
                        sequences[header] = ''.join(current_sequence)
                    
                    # Start a new sequence
                    header = line[1:] # Store header without the '>'
                    current_sequence = []
                elif header:
                    # Append sequence line
                    current_sequence.append(line)
            
            # Don't forget to save the very last sequence in the file
            if header:
                sequences[header] = ''.join(current_sequence)
                
    except:
        print(f"Error: The file '{file_path}' was not found.", file=sys.stderr)
        return None
        
    return sequences

def compare_sequences(dict1, dict2, file1_name, file2_name):
    """
    Compares two sequence dictionaries and prints a detailed report.
    """
    # Python's dictionary comparison is a powerful shortcut!
    # It checks for same keys and same values for each key.
    if dict1 == dict2:
        print("✅ The two FASTA files contain the exact same set of sequences and headers.")
        print(f"   - Found {len(dict1)} sequences in both files.")
        return True

    print("❌ The two FASTA files are different. Here's a summary:\n")
    
    # Get the set of headers from each dictionary
    headers1 = set(dict1.keys())
    headers2 = set(dict2.keys())

    # Find headers unique to each file
    unique_to_1 = headers1 - headers2
    unique_to_2 = headers2 - headers1
    
    # Find headers that are common to both files
    common_headers = headers1 & headers2

    # Check for differences in content for common headers
    differing_content = []
    for header in common_headers:
        if dict1[header] != dict2[header]:
            differing_content.append(header)
    
    # --- Generate Report ---
    is_different = False

    # Report on total sequence counts
    if len(dict1) != len(dict2):
        is_different = True
        print(f"Mismatched sequence counts:")
        print(f"  - '{file1_name}' has {len(dict1)} sequences.")
        print(f"  - '{file2_name}' has {len(dict2)} sequences.\n")

    # Report on headers unique to file 1
    if unique_to_1:
        is_different = True
        print(f"Headers found only in '{file1_name}' ({len(unique_to_1)}):")
        # Print first 5 for brevity
        for header in list(unique_to_1)[:5]:
            print(f"  - {header}")
        if len(unique_to_1) > 5:
            print(f"  - ... and {len(unique_to_1) - 5} more.")
        print()

    # Report on headers unique to file 2
    if unique_to_2:
        is_different = True
        print(f"Headers found only in '{file2_name}' ({len(unique_to_2)}):")
        # Print first 5 for brevity
        for header in list(unique_to_2)[:5]:
            print(f"  - {header}")
        if len(unique_to_2) > 5:
            print(f"  - ... and {len(unique_to_2) - 5} more.")
        print()

    # Report on sequences with same header but different content
    if differing_content:
        is_different = True
        print(f"Headers with different sequences ({len(differing_content)}):")
        # Print first 5 for brevity
        for header in differing_content[:5]:
            print(f"  - {header}")
        if len(differing_content) > 5:
            print(f"  - ... and {len(differing_content) - 5} more.")
        print()

    return not is_different

def eval_fasta_from_file_equal(res_data, ref_data):
    """Main function to run the script."""

    print(f"Comparing '{res_data}' and '{ref_data}'...\n")

    # Step 1: Parse both files into dictionaries
    seqs1 = parse_fasta(res_data)
    seqs2 = parse_fasta(ref_data)

    if seqs1 is None or seqs2 is None:
        return False
    
    # Step 2: Compare the dictionaries and report the results
    return compare_sequences(seqs1, seqs2, res_data, ref_data)

def _check_files_exist(file1_path: str, file2_path: str) -> bool:
    """辅助函数：检查两个文件是否存在，如果不存在则记录错误。"""
    if not os.path.exists(file1_path):
        print(f"File not found: {file1_path}")
        return False
    if not os.path.exists(file2_path):
        print(f"File not found: {file2_path}")
        return False
    return True


# --- TSV / Tabular Data Functions ---

def eval_tsv_file_sort(file1_path: str, file2_path: str, sep='\t', header='infer') -> bool:
    """
    评估两个TSV文件内容是否完全相同（包括行顺序）。

    使用 pandas 库进行健壮的比较，能正确处理数据类型、缺失值等。

    Args:
        file1_path (str): 第一个文件的路径。
        file2_path (str): 第二个文件的路径。
        sep (str): 字段分隔符，默认为制表符。
        header (int, 'infer', or None): 指定表头行。

    Returns:
        bool: 如果文件内容按顺序完全相同，则返回 True，否则返回 False。
    """
    if not _check_files_exist(file1_path, file2_path):
        return False

    try:
        import pandas as pd
    except:
        print("pandas is not installed. Please run 'pip install pandas'.")
        return False

    try:
        df1 = pd.read_csv(file1_path, sep=sep, header=header, keep_default_na=False)
        df2 = pd.read_csv(file2_path, sep=sep, header=header, keep_default_na=False)
        return df1.equals(df2)
    except Exception as e:
        print(f"An error occurred while comparing TSV files: {e}")
        return False


def eval_tsv_file_unsort(file1_path: str, file2_path: str, sep='\t', header='infer') -> bool:
    """
    评估两个TSV文件内容是否相同（忽略行顺序）。

    此函数会读取两个文件，按所有列对内容进行排序，然后进行比较。
    这确保了即使行顺序不同，只要内容集合相同，结果也为 True。

    Args:
        file1_path (str): 第一个文件的路径。
        file2_path (str): 第二个文件的路径。
        sep (str): 字段分隔符，默认为制表符。
        header (int, 'infer', or None): 指定表头行。

    Returns:
        bool: 如果文件内容集合相同（忽略顺序），则返回 True，否则返回 False。
    """
    if not _check_files_exist(file1_path, file2_path):
        return False

    try:
        import pandas as pd
    except:
        print("pandas is not installed. Please run 'pip install pandas'.")
        return False

    try:
        df1 = pd.read_csv(file1_path, sep=sep, header=header, keep_default_na=False)
        df2 = pd.read_csv(file2_path, sep=sep, header=header, keep_default_na=False)

        # 如果列数或行数不同，它们肯定不相等
        if df1.shape != df2.shape:
            return False
            
        # 如果没有列，则无法排序，直接比较
        if df1.shape[1] == 0:
            return df1.equals(df2)

        # 按所有列排序以创建规范化的表示
        df1_sorted = df1.sort_values(by=list(df1.columns)).reset_index(drop=True)
        df2_sorted = df2.sort_values(by=list(df2.columns)).reset_index(drop=True)

        return df1_sorted.equals(df2_sorted)
    except Exception as e:
        print(f"An error occurred while comparing unsorted TSV files: {e}")
        return False


# --- Bioinformatics File Functions ---

def eval_bed_file_equal(file1_path: str, file2_path: str) -> bool:
    """
    评估两个BED文件内容是否完全相同（包括行顺序）。

    BED文件是一种特殊的TSV格式，通常没有表头。
    此函数执行严格的、按顺序的比较。

    Args:
        file1_path (str): 第一个BED文件的路径。
        file2_path (str): 第二个BED文件的路径。

    Returns:
        bool: 如果文件内容按顺序完全相同，则返回 True，否则返回 False。
    """
    # BED文件是TSV的一种，可以直接调用 `eval_tsv_file_sort` 并指定没有表头
    return eval_tsv_file_sort(file1_path, file2_path, sep='\t', header=None)


def eval_vcf_file_equal(file1_path: str, file2_path: str) -> bool:
    """
    评估两个VCF文件内容是否完全相同。

    此函数会逐行比较，但对元信息行（以'##'开头）和数据行进行区分。
    它会忽略元信息行的顺序，但要求数据行（包括表头 '#CHROM'）的顺序完全一致。

    Args:
        file1_path (str): 第一个VCF文件的路径。
        file2_path (str): 第二个VCF文件的路径。

    Returns:
        bool: 如果文件的元信息集合相同且数据部分完全相同，则返回 True。
    """
    if not _check_files_exist(file1_path, file2_path):
        return False

    try:
        meta1, data1 = set(), []
        meta2, data2 = set(), []

        with open(file1_path, 'r') as f1:
            for line in f1:
                if line.startswith('##'):
                    meta1.add(line.strip())
                else:
                    data1.append(line.strip())

        with open(file2_path, 'r') as f2:
            for line in f2:
                if line.startswith('##'):
                    meta2.add(line.strip())
                else:
                    data2.append(line.strip())

        # 比较元信息（忽略顺序）和数据（保持顺序）
        return (meta1 == meta2) and (data1 == data2)
    except Exception as e:
        print(f"An error occurred while comparing VCF files: {e}")
        return False


# --- Complex Binary/Serialized Object Functions ---

def eval_rds_file_equal(file1_path: str, file2_path: str) -> bool:
    """
    评估两个RDS文件（R语言对象）在语义上是否相等。

    使用 rpy2 库加载R对象，并调用R的 `all.equal` 函数进行深度比较。
    这比二进制比较更可靠，因为它不关心序列化的微小差异。
    需要安装 R、rpy2 库 (`pip install rpy2`)。

    Args:
        file1_path (str): 第一个RDS文件的路径。
        file2_path (str): 第二个RDS文件的路径。

    Returns:
        bool: 如果R对象在语义上相等，则返回 True，否则返回 False。
    """
    if not _check_files_exist(file1_path, file2_path):
        return False

    try:
        from rpy2.robjects import r
        import rpy2.robjects as robjects
    except:
        print("rpy2 is not installed. Please run 'pip install rpy2'.")
        print("You also need a working R installation on your system.")
        return False

    try:
        # 加载R的 `readRDS` 和 `all.equal` 函数
        readRDS = r['readRDS']
        all_equal = r['all.equal']

        # 读取RDS文件为R对象
        obj1 = readRDS(file1_path)
        obj2 = readRDS(file2_path)

        # 使用 all.equal 进行比较
        result = all_equal(obj1, obj2)
        
        # R的 all.equal 在相等时返回 R 的 TRUE，否则返回描述差异的字符串
        return result[0] is robjects.vectors.BoolVector(True)[0]
    except Exception as e:
        print(f"An error occurred while comparing RDS files: {e}")
        return False


def eval_h5ad_file_equal(file1_path: str, file2_path: str, tolerance: float = 1e-6) -> bool:
    """
    评估两个H5AD文件（AnnData对象）在语义上是否相等。

    加载两个AnnData对象并比较它们的关键组件：
    - X 矩阵 (数据)
    - obs DataFrame (观测/细胞元数据)
    - var DataFrame (变量/基因元数据)
    对于浮点数数据，会使用指定的容差进行比较。
    需要安装 anndata 和其依赖 (`pip install anndata pandas`).

    Args:
        file1_path (str): 第一个 H5AD 文件的路径。
        file2_path (str): 第二个 H5AD 文件的路径。
        tolerance (float): 比较浮点数数据时的绝对和相对容差。

    Returns:
        bool: 如果 AnnData 对象的关键组件相等，则返回 True。
    """
    if not _check_files_exist(file1_path, file2_path):
        return False

    try:
        import anndata as ad
        import pandas as pd
        from scipy.sparse import spmatrix
    except:
        print("anndata or pandas is not installed. Please run 'pip install anndata pandas'.")
        return False

    try:
        adata1 = ad.read_h5ad(file1_path)
        adata2 = ad.read_h5ad(file2_path)

        # 1. 比较 .obs (细胞元数据)
        if not adata1.obs.equals(adata2.obs):
            print("AnnData comparison failed: .obs DataFrames are different.")
            return False

        # 2. 比较 .var (基因元数据)
        if not adata1.var.equals(adata2.var):
            print("AnnData comparison failed: .var DataFrames are different.")
            return False

        # 3. 比较 .X (主数据矩阵)
        x1, x2 = adata1.X, adata2.X
        
        # 检查是否一个是稀疏矩阵而另一个不是
        if isinstance(x1, spmatrix) != isinstance(x2, spmatrix):
            print("AnnData comparison failed: .X matrices have different sparsity.")
            return False
            
        if isinstance(x1, spmatrix): # 如果是稀疏矩阵
            # 比较非零元素的差异
            if (x1 - x2).nnz > 0:
                print("AnnData comparison failed: Sparse .X matrices are different.")
                return False
        else: # 如果是密集矩阵 (numpy array)
            if not np.allclose(x1, x2, rtol=tolerance, atol=tolerance):
                print(f"AnnData comparison failed: Dense .X matrices are different within tolerance {tolerance}.")
                return False

        # 注意：一个完整的比较可能还需要检查 .uns, .obsm, .varm, .layers 等
        # 这里只比较了最核心的三个组件，对于大多数用例已经足够。
        print("AnnData objects appear to be equal based on .X, .obs, and .var.")
        return True
    except Exception as e:
        print(f"An error occurred while comparing H5AD files: {e}")
        return False
    

import gzip
from typing import TextIO, Set, List

# --- 为了支持 VCF 和 VCF.GZ 的比较，我们重构一下 VCF 解析逻辑 ---

def _parse_vcf_stream(vcf_stream: TextIO) -> Tuple[Set[str], List[str]]:
    """
    从一个文本流（文件或解压后的流）中解析VCF内容。

    Args:
        vcf_stream: 一个可迭代的文本流，产生VCF文件的行。

    Returns:
        一个元组，包含：
        - 一个包含所有元信息行（##开头）的集合（无序）。
        - 一个包含表头和所有数据行的列表（有序）。
    """
    meta_lines = set()
    data_lines = []
    for line in vcf_stream:
        stripped_line = line.strip()
        if not stripped_line:  # 跳过空行
            continue
        if stripped_line.startswith('##'):
            meta_lines.add(stripped_line)
        else:
            data_lines.append(stripped_line)
    return meta_lines, data_lines


def eval_vcf_file_equal(file1_path: str, file2_path: str) -> bool:
    """
    (重构后) 评估两个VCF文件内容是否完全相同。

    此函数会逐行比较，但对元信息行（以'##'开头）和数据行进行区分。
    它会忽略元信息行的顺序，但要求数据行（包括表头 '#CHROM'）的顺序完全一致。

    Args:
        file1_path (str): 第一个VCF文件的路径。
        file2_path (str): 第二个VCF文件的路径。

    Returns:
        bool: 如果文件的元信息集合相同且数据部分完全相同，则返回 True。
    """
    if not _check_files_exist(file1_path, file2_path):
        return False

    try:
        with open(file1_path, 'r', encoding='utf-8') as f1:
            meta1, data1 = _parse_vcf_stream(f1)

        with open(file2_path, 'r', encoding='utf-8') as f2:
            meta2, data2 = _parse_vcf_stream(f2)

        # 比较元信息（忽略顺序）和数据（保持顺序）
        return (meta1 == meta2) and (data1 == data2)
    except Exception as e:
        print(f"An error occurred while comparing VCF files: {e}")
        return False

# --- 以下是您要求的 5 个新函数 ---

def eval_fastq_file_equal(file1_path: str, file2_path: str) -> bool:
    """
    评估两个 FASTQ 文件内容是否完全相同。

    FASTQ 文件的记录顺序非常重要，因此该函数执行严格的、逐字节的文本内容比较。
    它实际上是 `eval_str_from_file_equal` 的一个特定别名，以提高代码可读性。

    Args:
        file1_path (str): 第一个 FASTQ 文件的路径。
        file2_path (str): 第二个 FASTQ 文件的路径。

    Returns:
        bool: 如果文件内容完全相同，则返回 True，否则返回 False。
    """
    return eval_str_from_file_equal(file1_path, file2_path)


def eval_fq_gz_file_equal(file1_path: str, file2_path: str) -> bool:
    """
    评估两个 gzipped FASTQ (.fq.gz) 文件的解压后内容是否完全相同。

    此函数会动态解压两个文件并比较其文本内容，而不是比较压缩后的二进制文件。
    这确保了即使压缩级别或时间戳不同，只要原始 FASTQ 内容相同，结果也为 True。

    Args:
        file1_path (str): 第一个 .fq.gz 文件的路径。
        file2_path (str): 第二个 .fq.gz 文件的路径。

    Returns:
        bool: 如果解压后的文件内容完全相同，则返回 True，否则返回 False。
    """
    if not _check_files_exist(file1_path, file2_path):
        return False
    
    try:
        with gzip.open(file1_path, 'rt', encoding='utf-8') as f1:
            content1 = f1.read()
        with gzip.open(file2_path, 'rt', encoding='utf-8') as f2:
            content2 = f2.read()
        
        return content1 == content2
    # except gzip.BadGzipFile as e:
    #     print(f"Error: A file is not a valid gzip file. {e}")
    #     return False
    except Exception as e:
        print(f"An error occurred while comparing gzipped FASTQ files: {e}")
        return False


def eval_vcf_gz_file_equal(file1_path: str, file2_path: str) -> bool:
    """
    评估两个 gzipped VCF (.vcf.gz) 文件的解压后内容是否相同。

    该函数使用与 `eval_vcf_file_equal` 相同的逻辑：
    - 忽略元信息行（##开头）的顺序。
    - 要求数据行（包括#CHROM表头）的顺序严格一致。
    它通过动态解压文件流来实现这一点。

    Args:
        file1_path (str): 第一个 .vcf.gz 文件的路径。
        file2_path (str): 第二个 .vcf.gz 文件的路径。

    Returns:
        bool: 如果解压后文件的元信息集合相同且数据部分完全相同，则返回 True。
    """
    if not _check_files_exist(file1_path, file2_path):
        return False

    try:
        with gzip.open(file1_path, 'rt', encoding='utf-8') as f1:
            meta1, data1 = _parse_vcf_stream(f1)

        with gzip.open(file2_path, 'rt', encoding='utf-8') as f2:
            meta2, data2 = _parse_vcf_stream(f2)

        # 比较元信息（忽略顺序）和数据（保持顺序）
        return (meta1 == meta2) and (data1 == data2)
    # except gzip.BadGzipFile as e:
    #     print(f"Error: A file is not a valid gzip file. {e}")
    #     return False
    except Exception as e:
        print(f"An error occurred while comparing gzipped VCF files: {e}")
        return False


def eval_gtf_file_equal(file1_path: str, file2_path: str) -> bool:
    """
    评估两个 GTF 文件内容是否相同，忽略行顺序和注释。

    GTF (Gene Transfer Format) 文件是带有注释（以'#'开头）的制表符分隔文件。
    在功能上，行的顺序通常不重要。此函数使用 pandas 读取数据，忽略注释和行序，
    然后比较数据内容的集合是否相同。

    Args:
        file1_path (str): 第一个 GTF 文件的路径。
        file2_path (str): 第二个 GTF 文件的路径。

    Returns:
        bool: 如果两个文件包含相同的非注释行集合，则返回 True。
    """
    if not _check_files_exist(file1_path, file2_path):
        return False
        
    try:
        # 读取 GTF 文件，它是没有表头的 TSV，且注释以 '#' 开头
        df1 = pd.read_csv(file1_path, sep='\t', header=None, comment='#')
        df2 = pd.read_csv(file2_path, sep='\t', header=None, comment='#')

        # 如果行数或列数不同，则文件不同
        if df1.shape != df2.shape:
            return False
            
        # 如果文件为空，则它们相等
        if df1.empty:
            return True

        # 按所有列对数据进行排序，以忽略原始行顺序
        df1_sorted = df1.sort_values(by=list(df1.columns)).reset_index(drop=True)
        df2_sorted = df2.sort_values(by=list(df2.columns)).reset_index(drop=True)

        return df1_sorted.equals(df2_sorted)
    except Exception as e:
        print(f"An error occurred while comparing GTF files: {e}")
        return False


def eval_psl_file_equal(file1_path: str, file2_path: str) -> bool:
    """
    评估两个 PSL 文件内容是否相同，忽略行顺序和标准头部。

    PSL (PSL format) 文件是 BLAT 的输出格式，通常有5行文件头。
    此函数假定并跳过这5行，然后比较剩余的比对记录。由于比对记录的顺序
    通常不重要，函数会忽略行序进行比较。

    Args:
        file1_path (str): 第一个 PSL 文件的路径。
        file2_path (str): 第二个 PSL 文件的路径。

    Returns:
        bool: 如果两个文件包含相同的比对记录集合，则返回 True。
    """
    if not _check_files_exist(file1_path, file2_path):
        return False

    try:
        # PSL 文件是 TSV，没有列名，通常有5行头部需要跳过
        df1 = pd.read_csv(file1_path, sep='\t', header=None, skiprows=5)
        df2 = pd.read_csv(file2_path, sep='\t', header=None, skiprows=5)

        if df1.shape != df2.shape:
            return False

        if df1.empty:
            return True

        # 排序以忽略原始行顺序
        df1_sorted = df1.sort_values(by=list(df1.columns)).reset_index(drop=True)
        df2_sorted = df2.sort_values(by=list(df2.columns)).reset_index(drop=True)

        return df1_sorted.equals(df2_sorted)
    except Exception as e:
        print(f"An error occurred while comparing PSL files: {e}")
        return False