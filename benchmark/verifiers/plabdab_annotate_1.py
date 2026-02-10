import pandas as pd
from typing import Tuple, Dict


def plabdab_annotate_1(d) -> bool:
    is_perfect_match, metrics = compare_cdr_csv_files(
        ref_csv_path = 'benchmark/gold_results/plabdab_annotate_1.csv', 
        output_csv_path = f"{d}/plabdab_annotate_1.csv", 
    )
    return is_perfect_match

def compare_cdr_csv_files(ref_csv_path: str, output_csv_path: str) -> Tuple[bool, Dict]:
    """
    Compare two CSV files containing CDR sequence data.
    
    Columns expected:
    - "Accession": The Accession ID from the FASTA header.
    - "Numbered_Sequence": The sequence string with IMGT numbering gaps.
    - "CDR_Lengths": The formatted length string (e.g., 8_8_12).
    - "Chain_Type": The identified chain (either H or L).
    
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
        required_columns = ["Accession", "Numbered_Sequence", "CDR_Lengths", "Chain_Type"]
        
        for df, name in [(ref_df, "reference"), (output_df, "output")]:
            missing_cols = set(required_columns) - set(df.columns)
            if missing_cols:
                return False, {"error": f"{name} CSV missing columns: {missing_cols}"}
        
        # Set Accession as index for easy lookup
        ref_df = ref_df.set_index("Accession")
        output_df = output_df.set_index("Accession")
        
        # Get all unique Accession IDs
        all_accessions = set(ref_df.index).union(set(output_df.index))
        ref_accessions = set(ref_df.index)
        output_accessions = set(output_df.index)
        
        # Calculate matches
        true_positives = 0
        false_positives = 0
        false_negatives = 0
        
        matched_accessions = []
        mismatched_details = []
        
        # Check output sequences against reference (precision calculation)
        for accession in output_accessions:
            if accession in ref_accessions:
                # Get rows for this accession
                output_row = output_df.loc[accession]
                ref_row = ref_df.loc[accession]
                
                # Check if all three columns match
                columns_match = True
                mismatched_cols = []
                
                for col in ["Numbered_Sequence", "CDR_Lengths", "Chain_Type"]:
                    output_val = str(output_row[col]) if pd.notna(output_row[col]) else ""
                    ref_val = str(ref_row[col]) if pd.notna(ref_row[col]) else ""
                    
                    if output_val != ref_val:
                        columns_match = False
                        mismatched_cols.append(col)
                
                if columns_match:
                    true_positives += 1
                    matched_accessions.append(accession)
                else:
                    false_positives += 1
                    mismatched_details.append({
                        "accession": accession,
                        "reason": "column_mismatch",
                        "mismatched_columns": mismatched_cols,
                        "output_values": {col: str(output_row[col]) for col in required_columns[1:]},
                        "ref_values": {col: str(ref_row[col]) for col in required_columns[1:]}
                    })
            else:
                false_positives += 1
                mismatched_details.append({
                    "accession": accession,
                    "reason": "accession_not_in_reference"
                })
        
        # Check reference sequences not in output (recall calculation)
        for accession in ref_accessions:
            if accession not in output_accessions:
                false_negatives += 1
                mismatched_details.append({
                    "accession": accession,
                    "reason": "accession_missing_in_output"
                })
        
        # Calculate metrics
        precision = true_positives / (true_positives + false_positives) if (true_positives + false_positives) > 0 else 0
        recall = true_positives / (true_positives + false_negatives) if (true_positives + false_negatives) > 0 else 0
        
        if precision + recall == 0:
            f1_score = 0
        else:
            f1_score = 2 * precision * recall / (precision + recall)
        
        # Check if F1 score equals 1 (perfect match)
        is_perfect_match = abs(f1_score - 1.0) < 1e-10
        
        metrics = {
            "f1_score": f1_score,
            "precision": precision,
            "recall": recall,
            "true_positives": true_positives,
            "false_positives": false_positives,
            "false_negatives": false_negatives,
            "total_ref_entries": len(ref_df),
            "total_output_entries": len(output_df),
            "unique_ref_accessions": len(ref_accessions),
            "unique_output_accessions": len(output_accessions),
            "matched_accessions_count": len(matched_accessions),
            "mismatched_count": len(mismatched_details),
            "is_perfect_match": is_perfect_match,
            "matched_accessions_sample": matched_accessions[:20],  # Sample for readability
            "mismatched_details_sample": mismatched_details[:20]   # Sample for readability
        }
        
        return is_perfect_match, metrics
        
    except FileNotFoundError as e:
        return False, {"error": f"File not found: {str(e)}"}
    except pd.errors.EmptyDataError:
        return False, {"error": "CSV file is empty"}
    except Exception as e:
        return False, {"error": f"Error processing CSV files: {str(e)}"}
