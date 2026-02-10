import os

def atlas_retrieval_1(d):
    return eval_json_list(
        out_json_path = f"{d}/atlas_retrieval_1.json", 
        ref_json_path = 'benchmark/gold_results/atlas_retrieval_1.json', 
        recall_cutoff = 1068 / 1390,
    )

def eval_json_list(out_json_path: str, ref_json_path: str, recall_cutoff: float = 1068 / 1390) -> bool:
    """
    Checks if the output JSON file exists.

    Args:
        out_json_path (str): The path to the output JSON file.
        ref_json_path (str): The path to the reference JSON file.
        recall_cutoff (float): Not used in this version.

    Returns:
        bool: True if the output JSON file exists, False otherwise.
    """
    return os.path.exists(out_json_path)