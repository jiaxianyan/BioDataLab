import string
import random
import json
from typing import Dict, List, Tuple, Union
from pathlib import Path
import numpy as np


def atlas_retrieval_1(d):
    return eval_json_list(
        out_json_path = f"{d}/atlas_retrieval_1.json", 
        ref_json_path = 'benchmark/gold_results/atlas_retrieval_1.json', 
        recall_cutoff = 1068 / 1390,
    )

def eval_json_list(out_json_path: str, ref_json_path: str, recall_cutoff: float = 1068 / 1390) -> bool:
    """
    Reads two JSON files, extracts the lists, and compares them.

    Args:
        out_json_path (str): The path to the output JSON file.
        ref_json_path (str): The path to the reference JSON file.

    Returns:
        bool: True if the more than `recall_cutoff` percent of list are considered in the reference list, False otherwise.
    """
    try:
        with open(out_json_path, 'r') as f1:
            out_list = json.load(f1)

        with open(ref_json_path, 'r') as f2:
            ref_list = json.load(f2)

    except FileNotFoundError as e:
        print(f"Error: File not found. Details: {e}")
        return False
    except json.JSONDecodeError as e:
        print(f"Error: Could not decode JSON. Make sure files are valid. Details: {e}")
        return False
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return False

    # Ensure that the loaded data from both files are actually lists
    if not isinstance(out_list, list) or not isinstance(ref_list, list):
        print("Error: One or both of the JSON files do not contain a list at the top level.")
        return False

    
    try:
        out_list = sorted(out_list)
        ref_list = sorted(ref_list)
    except TypeError:
        print("Error: The lists contain non-sortable items (like dictionaries).")
        print("Cannot perform an order-insensitive comparison using the default sort method.")
        # For complex cases like lists of dictionaries, a more advanced comparison
        # would be needed. For this function, we'll consider them not the same.
        return False

    duplicated_set = set()
    recalled_count = 0
    for x in out_list:
        try:
            pdbid, chainid = x.split('_')
        except ValueError as e:
            print(f'Element {x}: {str(e)}')
        
        x = f'{pdbid.lower()}_{chainid}'
        if x not in duplicated_set:
            duplicated_set.add(x)
            if x in ref_list:
                recalled_count += 1
    # print(f'{recalled_count} / {len(ref_list)} = {recalled_count / len(ref_list)}')
    return (recalled_count / len(ref_list)) >= recall_cutoff