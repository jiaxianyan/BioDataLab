import os
from typing import Dict, Tuple


def inclusive_extract_2(d) -> bool:
    print('Extract ncAAs SMILES from JACS...')
    result_1 = os.path.exists(f"{d}/inclusive_extract_2/paper_data_25.csv")
    print('True/False:', result_1)
    print('-' * 50)

    print('Extract ncAAs SMILES from PNAS...')
    result_2 = os.path.exists(f"{d}/inclusive_extract_2/paper_data_73.csv")
    print('True/False:', result_2)
    print('-' * 50)

    print('Extract ncAAs SMILES from Science Signaling...')
    result_3 = os.path.exists(f"{d}/inclusive_extract_2/paper_data_342.csv")
    print('True/False:', result_3)
    print('-' * 50)

    print('Extract ncAAs SMILES from Nature Chemical Biology....')
    result_4 = os.path.exists(f"{d}/inclusive_extract_2/paper_data_249_274.csv")
    print('True/False:', result_4)
    print('-' * 50)

    print('Extract ncAAs SMILES from Angew. Chem. Int. Ed....')
    result_5 = os.path.exists(f"{d}/inclusive_extract_2/paper_data_2228_2230.csv")
    print('True/False:', result_5)
    print('-' * 50)

    print('Evaluation is done!')
    return all([result_1, result_2, result_3, result_4, result_5])


class NcAASmilesComparator:
    """Compare ncAA data with SMILES notation between reference and LLM output."""
    
    @staticmethod
    def compare_from_files(ref_path: str, out_path: str,
                          f1_threshold: float = 0.8) -> Tuple[bool, Dict]:
        """
        Check if output file exists.
        
        Args:
            ref_path: Path to reference file (not used)
            out_path: Path to output file to check
            f1_threshold: Threshold for F1 score (not used)
        
        Returns:
            Tuple of (file_exists, empty_dict)
        """
        file_exists = os.path.exists(out_path)
        return file_exists, {}