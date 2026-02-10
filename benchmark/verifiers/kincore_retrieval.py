import string
import random
import json
from typing import Dict, List, Tuple, Union
from pathlib import Path
import numpy as np
import pandas as pd
from Bio import SearchIO


def kincore_retrieval(d) -> bool:
    return eval_json_list(
        out_json_path = f'{d}/kincore_retrieval.json', 
        ref_json_path = 'benchmark/gold_results/kincore_retrieval.json',
        f1_threshold = 0.8,
    )

def validate_format(data: List, expected_type: type) -> bool:
    """Validate that data is a list with all elements of expected_type."""
    if not isinstance(data, list):
        return False
    return all(isinstance(item, expected_type) for item in data)


def eval_json_list(out_json_path: str, 
                    ref_json_path: str, 
                    f1_threshold: float = 0.8) -> Tuple[bool, dict]:
    """
    Evaluate DOI list comparison with detailed metrics.
    
    Args:
        out_json_path: Path to LLM output JSON file
        ref_json_path: Path to ground truth JSON file
        f1_threshold: Minimum F1 value to consider successful
    
    Returns:
        Tuple of (success_flag, metrics_dict)
    """
    try:
        # Load and validate data
        with open(out_json_path, 'r') as f:
            out_data = json.load(f)
        with open(ref_json_path, 'r') as f:
            ref_data = json.load(f)
        
        # Validate formats
        if not validate_format(out_data, str):
            return False, {"error": "Invalid output format: expected list of strings"}
        if not validate_format(ref_data, str):
            return False, {"error": "Invalid reference format: expected list of strings"}
        
        # Calculate metrics
        out_set = set(out_data)
        ref_set = set(ref_data)
        
        true_positives = out_set.intersection(ref_set)
        false_positives = out_set - ref_set
        false_negatives = ref_set - out_set
        
        tp_count = len(true_positives)
        fp_count = len(false_positives)
        fn_count = len(false_negatives)
        
        precision = tp_count / (tp_count + fp_count) if (tp_count + fp_count) > 0 else 0
        recall = tp_count / (tp_count + fn_count) if (tp_count + fn_count) > 0 else 0
        
        if precision + recall == 0:
            f1_score = 0
        else:
            f1_score = 2 * precision * recall / (precision + recall)
        
        metrics = {
            "precision": precision,
            "recall": recall,
            "f1_score": f1_score,
            "true_positives": list(true_positives),
            "false_positives": list(false_positives),
            "false_negatives": list(false_negatives),
            "true_positives_count": tp_count,
            "false_positives_count": fp_count,
            "false_negatives_count": fn_count,
            "out_count": len(out_data),
            "ref_count": len(ref_data),
            "f1_threshold": f1_threshold,
            "meets_threshold": f1_score >= f1_threshold
        }
        
        return metrics["meets_threshold"]
        
    except json.JSONDecodeError as e:
        return False
    except FileNotFoundError as e:
        return False
    except Exception as e:
        return False