import os
from typing import Tuple

def inclusive_retrieval(d):
    meets_threshold, metrics = eval_json_list(
        out_json_path = f"{d}/inclusive_retrieval.json", 
        ref_json_path = 'benchmark/gold_results/inclusive_retrieval.json', 
        recall_threshold = 0.8,
    )
    return meets_threshold

def eval_json_list(out_json_path: str, 
                   ref_json_path: str,
                   recall_threshold: float = 0.8) -> Tuple[bool, dict]:
    """
    Check if the generated file exists.
    
    Args:
        out_json_path: Path to LLM output JSON file
        ref_json_path: Path to ground truth JSON file
        recall_threshold: Minimum recall value to consider successful
    
    Returns:
        Tuple of (success_flag, metrics_dict)
    """
    file_exists = os.path.exists(out_json_path)
    
    metrics = {
        "file_exists": file_exists,
        "out_json_path": out_json_path
    }
    
    return file_exists, metrics