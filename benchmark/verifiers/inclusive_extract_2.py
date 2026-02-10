import pandas as pd
import json
from typing import Dict, List, Tuple, Optional, Set
import numpy as np
from rdkit import Chem


def inclusive_extract_2(d) -> bool:
    print('Extract ncAAs SMILES from JACS...')
    result_1, _ = NcAASmilesComparator.compare_from_files(ref_path='benchmark/gold_results/iNClusive/paper_data_25_ref.csv', out_path=f"{d}/inclusive_extract_2/paper_data_25.csv")
    print('True/False:', result_1)
    print('-' * 50)

    print('Extract ncAAs SMILES from PNAS...')
    result_2, _ = NcAASmilesComparator.compare_from_files(ref_path='benchmark/gold_results/iNClusive/paper_data_73_ref.csv', out_path=f"{d}/inclusive_extract_2/paper_data_73.csv")
    print('True/False:', result_2)
    print('-' * 50)

    print('Extract ncAAs SMILES from Science Signaling...')
    result_3, _ = NcAASmilesComparator.compare_from_files(ref_path='benchmark/gold_results/iNClusive/paper_data_342_ref.csv', out_path=f"{d}/inclusive_extract_2/paper_data_342.csv")
    print('True/False:', result_3)
    print('-' * 50)

    print('Extract ncAAs SMILES from Nature Chemical Biology....')
    result_4, _ = NcAASmilesComparator.compare_from_files(ref_path='benchmark/gold_results/iNClusive/paper_data_249_274_ref.csv', out_path=f"{d}/inclusive_extract_2/paper_data_249_274.csv")
    print('True/False:', result_4)
    print('-' * 50)

    print('Extract ncAAs SMILES from Angew. Chem. Int. Ed....')
    result_5, _ = NcAASmilesComparator.compare_from_files(ref_path='benchmark/gold_results/iNClusive/paper_data_2228_2230_ref.csv', out_path=f"{d}/inclusive_extract_2/paper_data_2228_2230.csv")
    print('True/False:', result_5)
    print('-' * 50)

    print('Evaluation is done!')
    return all([result_1, result_2, result_3, result_4, result_5])

class NcAASmilesComparator:
    """Compare ncAA data with SMILES notation between reference and LLM output."""
    
    REQUIRED_COLUMNS = [
        "ncAA abbreviation(s) used in the publication",
        "ncAA name, as mentioned in the publication",
        "ncAA SMILES notation"
    ]
    
    INDEX_COLUMN = "ncAA abbreviation(s) used in the publication"
    
    @staticmethod
    def validate_dataframe(df: pd.DataFrame) -> bool:
        """Validate dataframe has all required columns."""
        if not isinstance(df, pd.DataFrame):
            return False
        
        missing_columns = set(NcAASmilesComparator.REQUIRED_COLUMNS) - set(df.columns)
        return len(missing_columns) == 0
    
    @staticmethod
    def normalize_string(value: str) -> str:
        """Normalize string for comparison."""
        if pd.isna(value):
            return ""
        return str(value).strip().lower()
    
    @staticmethod
    def normalize_abbreviation(value: str) -> str:
        """
        Normalize ncAA abbreviation for comparison.
        Handle cases like 'not available', multiple abbreviations, etc.
        """
        if pd.isna(value):
            return ""
        
        value = str(value).strip().lower()
        if value == "not available":
            return value
        
        # Remove extra spaces and standardize formatting
        value = ' '.join(value.split())
        return value
    
    @staticmethod
    def normalize_name(value: str) -> str:
        """Normalize ncAA name for comparison."""
        if pd.isna(value):
            return ""
        
        value = str(value).strip().lower()
        if value == "not available":
            return value
        
        # Remove extra spaces and standardize formatting
        value = ' '.join(value.split())
        return value
    
    @staticmethod
    def get_canonical_smiles(smiles: str) -> Optional[str]:
        """
        Convert SMILES to canonical SMILES using RDKit.
        Returns None if SMILES is invalid or can't be parsed.
        """
        if pd.isna(smiles) or not isinstance(smiles, str):
            return None
        
        smiles = smiles.strip()
        if not smiles or smiles.lower() == "not available":
            return None
        
        try:
            # Try to parse SMILES
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                # Try with sanitization turned off first
                mol = Chem.MolFromSmiles(smiles, sanitize=False)
                if mol is None:
                    return None
                
                # Try to sanitize
                try:
                    Chem.SanitizeMol(mol)
                except:
                    # If sanitization fails, use as-is but canonicalize
                    pass
            
            # Generate canonical SMILES
            canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
            return canonical_smiles
            
        except Exception as e:
            # If RDKit fails, fall back to string normalization
            print(f"RDKit SMILES parsing error for '{smiles}': {e}")
            return smiles.strip()
    
    @staticmethod
    def smiles_match(smiles1: str, smiles2: str) -> bool:
        """
        Compare two SMILES strings using canonical forms.
        Returns True if they represent the same molecule.
        """
        # Handle special cases
        if pd.isna(smiles1) or pd.isna(smiles2):
            return False
        
        s1 = str(smiles1).strip()
        s2 = str(smiles2).strip()
        
        # Check if both are 'not available'
        if s1.lower() == "not available" and s2.lower() == "not available":
            return True
        
        # If one is 'not available' and the other isn't
        if (s1.lower() == "not available") != (s2.lower() == "not available"):
            return False
        
        # Get canonical SMILES
        canonical1 = NcAASmilesComparator.get_canonical_smiles(s1)
        canonical2 = NcAASmilesComparator.get_canonical_smiles(s2)
        
        # If canonicalization failed, fall back to string comparison
        if canonical1 is None and canonical2 is None:
            return s1 == s2
        elif canonical1 is None or canonical2 is None:
            return False
        
        # Compare canonical SMILES
        return canonical1 == canonical2
    
    @staticmethod
    def rows_match(row1: pd.Series, row2: pd.Series) -> bool:
        """Check if two rows match based on all columns."""
        # Check abbreviation
        abbrev1 = NcAASmilesComparator.normalize_abbreviation(
            row1["ncAA abbreviation(s) used in the publication"]
        ) if "ncAA abbreviation(s) used in the publication" in row1 else ""
        
        abbrev2 = NcAASmilesComparator.normalize_abbreviation(
            row2["ncAA abbreviation(s) used in the publication"]
        ) if "ncAA abbreviation(s) used in the publication" in row2 else ""
        
        if abbrev1 != abbrev2:
            return False
        
        # Check name
        name1 = NcAASmilesComparator.normalize_name(
            row1["ncAA name, as mentioned in the publication"]
        ) if "ncAA name, as mentioned in the publication" in row1 else ""
        
        name2 = NcAASmilesComparator.normalize_name(
            row2["ncAA name, as mentioned in the publication"]
        ) if "ncAA name, as mentioned in the publication" in row2 else ""
        
        if name1 != name2:
            return False
        
        # Check SMILES
        smiles1 = row1["ncAA SMILES notation"] if "ncAA SMILES notation" in row1 else None
        smiles2 = row2["ncAA SMILES notation"] if "ncAA SMILES notation" in row2 else None
        
        if not NcAASmilesComparator.smiles_match(smiles1, smiles2):
            return False
        
        return True
    
    @staticmethod
    def create_index_key(row: pd.Series) -> str:
        """Create unique key from index column."""
        if NcAASmilesComparator.INDEX_COLUMN in row:
            value = NcAASmilesComparator.normalize_abbreviation(
                row[NcAASmilesComparator.INDEX_COLUMN]
            )
            return value
        return ""
    
    @staticmethod
    def compare_dataframes(ref_df: pd.DataFrame, out_df: pd.DataFrame,
                          f1_threshold: float = 0.8) -> Tuple[bool, Dict]:
        """
        Compare reference and output dataframes.
        
        Args:
            ref_df: Reference dataframe
            out_df: Output dataframe
            f1_threshold: Threshold for F1 score to consider match successful
        
        Returns:
            Tuple of (match_success, metrics_dict)
        """
        # Validate dataframes
        if not NcAASmilesComparator.validate_dataframe(ref_df):
            return False, {"error": "Reference dataframe missing required columns"}
        if not NcAASmilesComparator.validate_dataframe(out_df):
            return False, {"error": "Output dataframe missing required columns"}
        
        # Create index dictionaries
        ref_dict = {}
        for idx, row in ref_df.iterrows():
            key = NcAASmilesComparator.create_index_key(row)
            if key:  # Only add non-empty keys
                ref_dict[key] = row
        
        out_dict = {}
        for idx, row in out_df.iterrows():
            key = NcAASmilesComparator.create_index_key(row)
            if key:  # Only add non-empty keys
                out_dict[key] = row
        
        # Calculate matches
        true_positives = 0
        false_positives = 0
        false_negatives = 0
        
        matched_rows = []
        mismatched_rows = []
        
        # Check output rows against reference
        for out_key, out_row in out_dict.items():
            if out_key in ref_dict:
                if NcAASmilesComparator.rows_match(out_row, ref_dict[out_key]):
                    true_positives += 1
                    matched_rows.append({
                        "key": out_key,
                        "out_row": out_row.to_dict(),
                        "ref_row": ref_dict[out_key].to_dict()
                    })
                else:
                    false_positives += 1
                    mismatched_rows.append({
                        "key": out_key,
                        "reason": "Row content mismatch",
                        "out_row": out_row.to_dict(),
                        "ref_row": ref_dict[out_key].to_dict()
                    })
            else:
                false_positives += 1
                mismatched_rows.append({
                    "key": out_key,
                    "reason": "Key not found in reference",
                    "out_row": out_row.to_dict()
                })
        
        # Check reference rows not in output
        for ref_key in ref_dict:
            if ref_key not in out_dict:
                false_negatives += 1
                mismatched_rows.append({
                    "key": ref_key,
                    "reason": "Key missing in output",
                    "ref_row": ref_dict[ref_key].to_dict()
                })
        
        # Calculate metrics
        precision = (true_positives / (true_positives + false_positives) 
                    if (true_positives + false_positives) > 0 else 0)
        recall = (true_positives / (true_positives + false_negatives) 
                 if (true_positives + false_negatives) > 0 else 0)
        
        if precision + recall == 0:
            f1_score = 0
        else:
            f1_score = 2 * precision * recall / (precision + recall)
        
        # Analyze SMILES-specific mismatches
        smiles_mismatches = []
        for mismatch in mismatched_rows:
            if "out_row" in mismatch and "ref_row" in mismatch:
                out_smiles = mismatch["out_row"].get("ncAA SMILES notation")
                ref_smiles = mismatch["ref_row"].get("ncAA SMILES notation")
                
                canonical_out = NcAASmilesComparator.get_canonical_smiles(out_smiles)
                canonical_ref = NcAASmilesComparator.get_canonical_smiles(ref_smiles)
                
                if canonical_out != canonical_ref:
                    smiles_mismatches.append({
                        "key": mismatch["key"],
                        "out_smiles": out_smiles,
                        "ref_smiles": ref_smiles,
                        "canonical_out": canonical_out,
                        "canonical_ref": canonical_ref
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
            "total_unique_keys_ref": len(ref_dict),
            "total_unique_keys_out": len(out_dict),
            "matched_rows_count": len(matched_rows),
            "mismatched_rows_count": len(mismatched_rows),
            "smiles_mismatches_count": len(smiles_mismatches),
            "f1_threshold": f1_threshold,
            "meets_threshold": f1_score >= f1_threshold,
            "matched_rows": matched_rows[:10],  # Limit for readability
            "smiles_mismatches": smiles_mismatches[:10]  # Limit for readability
        }
        
        return metrics["meets_threshold"], metrics
    
    @staticmethod
    def compare_from_files(ref_path: str, out_path: str,
                          f1_threshold: float = 0.8) -> Tuple[bool, Dict]:
        """
        Compare reference and output files.
        
        Args:
            ref_path: Path to reference file (CSV, Excel, or JSON)
            out_path: Path to output file (CSV, Excel, or JSON)
            f1_threshold: Threshold for F1 score to consider match successful
        
        Returns:
            Tuple of (match_success, metrics_dict)
        """
        try:
            # Load data based on file extension
            def load_file(filepath: str) -> pd.DataFrame:
                if filepath.endswith('.csv'):
                    return pd.read_csv(filepath, sep = ';')
                elif filepath.endswith(('.xlsx', '.xls')):
                    return pd.read_excel(filepath)
                elif filepath.endswith('.json'):
                    with open(filepath, 'r') as f:
                        data = json.load(f)
                    return pd.DataFrame(data)
                else:
                    # Assume CSV
                    return filepath
            
            ref_df = load_file(ref_path)
            out_df = load_file(out_path)
            
            return NcAASmilesComparator.compare_dataframes(ref_df, out_df, f1_threshold)
            
        except FileNotFoundError as e:
            return False, {"error": f"File not found: {str(e)}"}
        except json.JSONDecodeError as e:
            return False, {"error": f"JSON parsing error: {str(e)}"}
        except pd.errors.EmptyDataError:
            return False, {"error": "File is empty"}
        except Exception as e:
            return False, {"error": f"Unexpected error: {str(e)}"}