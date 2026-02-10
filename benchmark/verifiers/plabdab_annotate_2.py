import pandas as pd
from typing import Tuple, Dict, List, Set


def plabdab_annotate_2(d) -> bool:
    is_perfect_match, metrics = compare_antibody_csv_with_multiple_targets(
        ref_csv_path = 'benchmark/gold_results/plabdab_annotate_2.csv', 
        output_csv_path = f"{d}/plabdab_annotate_2.csv", 
    )
    return is_perfect_match

def compare_antibody_csv_with_multiple_targets(ref_csv_path: str, output_csv_path: str) -> Tuple[bool, Dict]:
    """
    Compare two CSV files containing antibody target information.
    Reference CSV may contain multiple correct targets separated by commas.
    
    Columns expected:
    - "ID": The Accession ID from the FASTA header.
    - "target_mention": Specific antigen being antibody-targeted.
        In reference: comma-separated list of correct targets
        In output: single target (matches if it matches any reference target)
    
    Args:
        ref_csv_path: Path to reference CSV file
        output_csv_path: Path to output CSV file
    
    Returns:
        Tuple of (is_match, metrics_dict)
        is_match: True if F1 score == 1 (perfect match), False otherwise
        metrics_dict: Dictionary containing detailed comparison metrics
    """
    try:
        # Load CSV files
        ref_df = pd.read_csv(ref_csv_path)
        output_df = pd.read_csv(output_csv_path)
        
        # Validate required columns exist
        required_columns = ["ID", "target_mention"]
        
        for df, name in [(ref_df, "reference"), (output_df, "output")]:
            missing_cols = set(required_columns) - set(df.columns)
            if missing_cols:
                return False, {"error": f"{name} CSV missing columns: {missing_cols}"}
        
        # Clean and prepare data
        # Remove rows with NaN in ID column
        ref_df = ref_df.dropna(subset=["ID"])
        output_df = output_df.dropna(subset=["ID"])
        
        # Convert ID to string for consistent comparison
        ref_df["ID"] = ref_df["ID"].astype(str).str.strip()
        output_df["ID"] = output_df["ID"].astype(str).str.strip()
        
        # Helper function to normalize and parse targets
        def parse_targets(target_str: str) -> Set[str]:
            """Parse comma-separated targets into a set of normalized strings."""
            if pd.isna(target_str):
                return set()
            
            targets = set()
            # Split by comma, handle potential spaces
            for target in str(target_str).split(','):
                target = target.strip().lower()
                if target:  # Only add non-empty targets
                    targets.add(target)
            return targets
        
        def normalize_target(target_str: str) -> str:
            """Normalize single target string."""
            if pd.isna(target_str):
                return ""
            return str(target_str).strip().lower()
        
        # Parse reference targets into sets
        ref_df["target_set"] = ref_df["target_mention"].apply(parse_targets)
        output_df["target_normalized"] = output_df["target_mention"].apply(normalize_target)
        
        # Create dictionaries for lookup
        ref_dict = ref_df.set_index("ID")[["target_set"]].to_dict(orient="index")
        output_dict = output_df.set_index("ID")[["target_normalized"]].to_dict(orient="index")
        
        # Get all unique IDs
        all_ids = set(ref_dict.keys()).union(set(output_dict.keys()))
        ref_ids = set(ref_dict.keys())
        output_ids = set(output_dict.keys())
        
        # Calculate matches
        true_positives = 0
        false_positives = 0
        false_negatives = 0
        
        matched_ids = []
        mismatched_details = []
        
        # Check output entries against reference (precision calculation)
        for id_val in output_ids:
            if id_val in ref_ids:
                # Get targets
                output_target = output_dict[id_val].get("target_normalized", "")
                ref_target_set = ref_dict[id_val].get("target_set", set())
                
                # Check if output target matches any in reference set
                if output_target in ref_target_set:
                    true_positives += 1
                    matched_ids.append(id_val)
                else:
                    false_positives += 1
                    mismatched_details.append({
                        "id": id_val,
                        "reason": "target_mismatch",
                        "output_target": output_target,
                        "ref_targets": list(ref_target_set),
                        "output_matches_ref": False
                    })
            else:
                false_positives += 1
                mismatched_details.append({
                    "id": id_val,
                    "reason": "id_not_in_reference"
                })
        
        # Check reference entries not in output (recall calculation)
        for id_val in ref_ids:
            if id_val not in output_ids:
                false_negatives += 1
                mismatched_details.append({
                    "id": id_val,
                    "reason": "id_missing_in_output",
                    "ref_targets": list(ref_dict[id_val].get("target_set", set()))
                })
        
        # Calculate metrics
        precision = true_positives / (true_positives + false_positives) if (true_positives + false_positives) > 0 else 0
        recall = true_positives / (true_positives + false_negatives) if (true_positives + false_negatives) > 0 else 0
        
        if precision + recall == 0:
            f1_score = 0
        else:
            f1_score = 2 * precision * recall / (precision + recall)
        
        # Check if F1 score equals 1 (perfect match)
        # Using tolerance for floating point comparison
        is_perfect_match = abs(f1_score - 1.0) < 1e-10
        
        # Calculate additional statistics
        total_ref_multiple_targets = sum(1 for targets in ref_df["target_set"] if len(targets) > 1)
        
        metrics = {
            "f1_score": f1_score,
            "precision": precision,
            "recall": recall,
            "true_positives": true_positives,
            "false_positives": false_positives,
            "false_negatives": false_negatives,
            "total_ref_entries": len(ref_df),
            "total_output_entries": len(output_df),
            "unique_ref_ids": len(ref_ids),
            "unique_output_ids": len(output_ids),
            "matched_ids_count": len(matched_ids),
            "mismatched_count": len(mismatched_details),
            "ref_entries_with_multiple_targets": total_ref_multiple_targets,
            "is_perfect_match": is_perfect_match,
            "matched_ids_sample": matched_ids[:20] if len(matched_ids) > 20 else matched_ids,
            "mismatched_details_sample": mismatched_details[:20] if len(mismatched_details) > 20 else mismatched_details
        }
        
        return is_perfect_match, metrics
        
    except FileNotFoundError as e:
        return False, {"error": f"File not found: {str(e)}"}
    except pd.errors.EmptyDataError:
        return False, {"error": "CSV file is empty"}
    except KeyError as e:
        return False, {"error": f"Missing expected column: {str(e)}"}
    except Exception as e:
        return False, {"error": f"Error processing CSV files: {str(e)}"}