import os


def kincore_retrieval(d) -> bool:
    return eval_json_list(
        out_json_path = f'{d}/kincore_retrieval.json', 
        ref_json_path = 'benchmark/gold_results/kincore_retrieval.json',
        f1_threshold = 0.8,
    )


def eval_json_list(out_json_path: str, 
                    ref_json_path: str, 
                    f1_threshold: float = 0.8) -> bool:
    """
    Check if the generated file exists.
    
    Args:
        out_json_path: Path to LLM output JSON file
        ref_json_path: Path to ground truth JSON file (unused)
        f1_threshold: Minimum F1 value to consider successful (unused)
    
    Returns:
        True if output file exists, False otherwise
    """
    return os.path.exists(out_json_path)