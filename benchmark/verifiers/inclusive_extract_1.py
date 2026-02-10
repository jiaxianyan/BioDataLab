import pandas as pd
import json
from typing import Dict, List, Tuple, Optional, Set, Union
import numpy as np


def inclusive_extract_1(d) -> bool:
    print('Retrieve ncAAs metadata from Journal of Translational Medicine...')
    result_1, _ = TableDataComparator.compare_tables_from_files('benchmark/gold_results/iNClusive/paper_data_2126_2127_ref.csv', f"{d}/inclusive_extract_1/paper_data_2126_2127.csv")
    print('True/False:', result_1)
    print('-' * 50)

    print('Retrieve ncAAs metadata from ACS Chemical Biology....')
    result_2, _ = TableDataComparator.compare_tables_from_files('benchmark/gold_results/iNClusive/paper_data_1748_ref.csv', f"{d}/inclusive_extract_1/paper_data_1748.csv")
    print('True/False:', result_2)
    print('-' * 50)

    print('Retrieve ncAAs metadata from PNAS....')
    result_3, _ = TableDataComparator.compare_tables_from_files('benchmark/gold_results/iNClusive/paper_data_947_ref.csv', f"{d}/inclusive_extract_1/paper_data_947.csv")
    print('True/False:', result_3)
    print('-' * 50)

    print('Retrieve ncAAs metadata from Nature Chemical Biology and its SI....')
    result_4, _ = TableDataComparator.compare_tables_from_files('benchmark/gold_results/iNClusive/paper_data_966_1011_ref.csv', f"{d}/inclusive_extract_1/paper_data_966_1011.csv")
    print('True/False:', result_4)
    print('-' * 50)

    print('Retrieve ncAAs metadata from Science....')
    result_5, _ = TableDataComparator.compare_tables_from_files('benchmark/gold_results/iNClusive/paper_data_17_21_ref.csv', f"{d}/inclusive_extract_1/paper_data_17_21.csv")
    print('True/False:', result_5)
    print('-' * 50)

    print('Evaluation is done!')
    return all([result_1, result_2, result_3, result_4, result_5])

class TableDataComparator:
    """Compare table data between reference and LLM output."""
    
    # Required columns for validation
    REQUIRED_COLUMNS = [
        "ncAA abbreviation(s) used in the publication",
        "ncAA name, as mentioned in the publication",
        "aaRS ID (abbr. organism, abbr. natural substrate, RS, mutations if any)",
        "aaRS origin organism full name",
        "tRNA ID (abbr. organism, tRNA, natural AA transported, anticodon)",
        "tRNA organism",
        "Tested in (protein)",
        "Tested in (protein position)",
        "Tested in (organism/in vitro)"
    ]
    
    # Index columns for matching rows
    INDEX_COLUMNS = [
        "ncAA abbreviation(s) used in the publication",
        "aaRS ID (abbr. organism, abbr. natural substrate, RS, mutations if any)",
        "tRNA ID (abbr. organism, tRNA, natural AA transported, anticodon)",
        "Tested in (protein)",
        "Tested in (protein position)"
    ]
    
    @staticmethod
    def validate_dataframe(df: pd.DataFrame) -> bool:
        """Validate dataframe has all required columns."""
        if not isinstance(df, pd.DataFrame):
            return False
        
        # Check if all required columns are present
        missing_columns = set(TableDataComparator.REQUIRED_COLUMNS) - set(df.columns)
        return len(missing_columns) == 0
    
    @staticmethod
    def normalize_string(value: str) -> str:
        """Normalize string for comparison."""
        if pd.isna(value):
            return ""
        return str(value).strip().lower()
    
    @staticmethod
    def normalize_position(value: str) -> str:
        """Normalize protein position values for comparison."""
        if pd.isna(value):
            return ""
        
        value = str(value).strip()
        if value == "not available":
            return value
        
        # Handle position formats
        if ";" in value:
            # Sort positions separated by semicolon
            positions = [p.strip() for p in value.split(";")]
            try:
                # Try to sort numerically
                positions.sort(key=lambda x: int(''.join(filter(str.isdigit, x)) or '0'))
            except:
                # If numeric parsing fails, sort alphabetically
                positions.sort()
            return ";".join(positions)
        return value
    
    @staticmethod
    def rows_match(row1: pd.Series, row2: pd.Series) -> bool:
        """Check if two rows match based on all columns."""
        for col in TableDataComparator.REQUIRED_COLUMNS:
            val1 = TableDataComparator.normalize_string(row1[col]) if col in row1 else ""
            val2 = TableDataComparator.normalize_string(row2[col]) if col in row2 else ""
            
            # Special handling for protein positions
            if col == "Tested in (protein position)":
                val1 = TableDataComparator.normalize_position(row1[col]) if col in row1 else ""
                val2 = TableDataComparator.normalize_position(row2[col]) if col in row2 else ""
            
            if val1 != val2:
                return False
        return True
    
    @staticmethod
    def create_row_index(row: pd.Series) -> str:
        """Create unique index string from index columns."""
        index_parts = []
        for col in TableDataComparator.INDEX_COLUMNS:
            if col in row:
                value = TableDataComparator.normalize_string(row[col])
                # Special handling for protein positions
                if col == "Tested in (protein position)":
                    value = TableDataComparator.normalize_position(row[col])
                index_parts.append(value)
            else:
                index_parts.append("")
        return "||".join(index_parts)
    
    @staticmethod
    def compare_tables_from_files(ref_path: str, out_path: str, 
                                 f1_threshold: float = 0.8) -> Tuple[bool, Dict]:
        """
        Compare reference and output tables.
        
        Args:
            ref_path: Path to reference CSV/Excel file
            out_path: Path to output CSV/Excel file
            f1_threshold: Threshold for F1 score to consider match successful
        
        Returns:
            Tuple of (match_success, metrics_dict)
        """
        try:
            # Load data based on file extension
            if ref_path.endswith('.csv'):
                ref_df = pd.read_csv(ref_path, sep = ';')
            elif ref_path.endswith(('.xlsx', '.xls')):
                ref_df = pd.read_excel(ref_path)
            else:
                # Assume pd.DataFrame
                ref_df = ref_path
            
            if out_path.endswith('.csv'):
                out_df = pd.read_csv(out_path, sep = ';')
            elif out_path.endswith(('.xlsx', '.xls')):
                out_df = pd.read_excel(out_path)
            else:
                # Assume pd.DataFrame
                out_df = out_path
            
            # Validate dataframes
            if not TableDataComparator.validate_dataframe(ref_df):
                return False, {"error": "Reference dataframe missing required columns"}
            if not TableDataComparator.validate_dataframe(out_df):
                return False, {"error": "Output dataframe missing required columns"}
            
            # Create index dictionaries for both dataframes
            ref_indices = {}
            for idx, row in ref_df.iterrows():
                row_index = TableDataComparator.create_row_index(row)
                ref_indices[row_index] = row
            
            out_indices = {}
            for idx, row in out_df.iterrows():
                row_index = TableDataComparator.create_row_index(row)
                out_indices[row_index] = row
            
            # Calculate matches
            true_positives = 0
            false_positives = 0
            false_negatives = 0
            
            # Check output rows against reference
            for out_index, out_row in out_indices.items():
                if out_index in ref_indices:
                    if TableDataComparator.rows_match(out_row, ref_indices[out_index]):
                        true_positives += 1
                    else:
                        false_positives += 1
                else:
                    false_positives += 1
            
            # Check reference rows not in output
            for ref_index in ref_indices:
                if ref_index not in out_indices:
                    false_negatives += 1
            
            # Calculate metrics
            precision = (true_positives / (true_positives + false_positives) 
                        if (true_positives + false_positives) > 0 else 0)
            recall = (true_positives / (true_positives + false_negatives) 
                     if (true_positives + false_negatives) > 0 else 0)
            
            if precision + recall == 0:
                f1_score = 0
            else:
                f1_score = 2 * precision * recall / (precision + recall)
            
            # Detailed match information
            matched_rows = []
            mismatched_rows = []
            
            for out_index, out_row in out_indices.items():
                if out_index in ref_indices:
                    if TableDataComparator.rows_match(out_row, ref_indices[out_index]):
                        matched_rows.append({
                            "index": out_index,
                            "row_data": out_row.to_dict()
                        })
                    else:
                        mismatched_rows.append({
                            "index": out_index,
                            "out_data": out_row.to_dict(),
                            "ref_data": ref_indices[out_index].to_dict()
                        })
            
            metrics = {
                "f1_score": f1_score,
                "precision": precision,
                "recall": recall,
                "true_positives": true_positives,
                "false_positives": false_positives,
                "false_negatives": false_negatives,
                "total_ref_rows": len(ref_df),
                "total_out_rows": len(out_df),
                "matched_rows_count": len(matched_rows),
                "mismatched_rows_count": len(mismatched_rows),
                "f1_threshold": f1_threshold,
                "meets_threshold": f1_score >= f1_threshold,
                "matched_rows": matched_rows,
                "mismatched_rows": mismatched_rows
            }
            
            return metrics["meets_threshold"], metrics
            
        except FileNotFoundError as e:
            return False, {"error": f"File not found: {str(e)}"}
        except pd.errors.EmptyDataError:
            return False, {"error": "File is empty"}
        except Exception as e:
            return False, {"error": f"Unexpected error: {str(e)}"}