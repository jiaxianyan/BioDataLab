import random
import json
from typing import Dict, List, Tuple, Union


def atlas_retrieval_2(d) -> bool:
    match, results = compare_json_files_with_details(
        out_json_path = f'{d}/atlas_retrieval_2.json', 
        ref_json_path = 'benchmark/gold_results/atlas_retrieval_2.json', 
    )
    return match


def compare_json_files_with_details(out_json_path: str, ref_json_path: str) -> Tuple[bool, Dict]:
    """
    Compare JSON files and return detailed results including F1 score components.
    Returns (F1_is_one, results_dict) where results_dict contains detailed metrics.
    """
    try:
        # Load JSON files
        with open(out_json_path, 'r') as f:
            out_data = json.load(f)
        with open(ref_json_path, 'r') as f:
            ref_data = json.load(f)
        
        # Validate out_data format
        if not isinstance(out_data, dict):
            return False, {"error": "Invalid output format: not a dictionary"}
        
        for key, value in out_data.items():
            if not isinstance(key, str) or not isinstance(value, dict):
                return False, {"error": "Invalid output format: keys must be strings and values must be dictionaries"}
            for inner_key, inner_value in value.items():
                if not isinstance(inner_key, str) or not isinstance(inner_value, str):
                    return False, {"error": "Invalid output format: inner keys must be strings and inner values must be strings"}
        
        # Calculate matches
        out_matches = {}
        ref_matches = {}
        
        # Out perspective
        for out_key, out_value in out_data.items():
            out_matches[out_key] = {}
            if out_key in ref_data:
                for inner_key, predicted_value in out_value.items():
                    if (inner_key in ref_data[out_key] and 
                        predicted_value in ref_data[out_key][inner_key]):
                        out_matches[out_key][inner_key] = True
                    else:
                        out_matches[out_key][inner_key] = False
            else:
                for inner_key in out_value.keys():
                    out_matches[out_key][inner_key] = False
        
        # Ref perspective
        for ref_key, ref_value in ref_data.items():
            ref_matches[ref_key] = {}
            if ref_key in out_data:
                for inner_key, true_values in ref_value.items():
                    if (inner_key in out_data[ref_key] and 
                        out_data[ref_key][inner_key] in true_values):
                        ref_matches[ref_key][inner_key] = True
                    else:
                        ref_matches[ref_key][inner_key] = False
            else:
                for inner_key in ref_value.keys():
                    ref_matches[ref_key][inner_key] = False
        
        # Calculate metrics
        total_out = sum(len(v) for v in out_data.values())
        total_ref = sum(len(v) for v in ref_data.values())
        
        tp_out = sum(sum(1 for v in inner.values() if v) for inner in out_matches.values())
        tp_ref = sum(sum(1 for v in inner.values() if v) for inner in ref_matches.values())
        
        precision = tp_out / total_out if total_out > 0 else 0
        recall = tp_ref / total_ref if total_ref > 0 else 0
        
        if precision + recall == 0:
            f1_score = 0
        else:
            f1_score = 2 * precision * recall / (precision + recall)
        
        results = {
            "f1_score": f1_score,
            "precision": precision,
            "recall": recall,
            "true_positives_out": tp_out,
            "true_positives_ref": tp_ref,
            "total_out_keys": total_out,
            "total_ref_keys": total_ref,
            "out_matches": out_matches,
            "ref_matches": ref_matches,
            "f1_is_one": abs(f1_score - 1.0) < 1e-10
        }
        
        return results["f1_is_one"], results
        
    except Exception as e:
        return False, {"error": str(e)}