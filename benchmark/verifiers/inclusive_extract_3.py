import json
from typing import Dict, Any


def inclusive_extract_3(d) -> bool:
    return compare_json_files(
        ref_json_path = 'benchmark/gold_results/iNClusive/inclusive_extract_3.json', 
        out_json_path = f"{d}/inclusive_extract_3.json",
    )

def compare_json_files(ref_json_path: str, out_json_path: str) -> bool:
    """
    Compare if two JSON files contain dict data.
    
    Args:
        ref_json_path: Path to reference JSON file
        out_json_path: Path to LLM output JSON file
    
    Returns:
        True if both files have identical content, False otherwise
    """
    try:
        # Load JSON files
        with open(ref_json_path, 'r') as f:
            ref_data = json.load(f)
        
        with open(out_json_path, 'r') as f:
            out_data = json.load(f)
        
        # Check if both are dictionaries
        if not isinstance(ref_data, dict) or not isinstance(out_data, dict):
            return False
        
        # Check if both have exactly the same keys
        if set(ref_data.keys()) != set(out_data.keys()):
            return False
        
        # Check if all values match exactly
        for key, ref_value in ref_data.items():
            out_value = out_data.get(key)
            
            # Special handling for string comparison (exact match required)
            if isinstance(ref_value, str) and isinstance(out_value, str):
                if ref_value != out_value:
                    return False
            else:
                # For non-string values or mixed types, use standard equality
                if ref_value != out_value:
                    return False
        
        return True
        
    except (json.JSONDecodeError, FileNotFoundError):
        # Handle JSON parsing errors or missing files
        return False
    except Exception:
        # Catch any other unexpected errors
        return False