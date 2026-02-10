import os
from typing import Tuple, Dict


def atlas_retrieval_2(d) -> bool:
    return check_file_exists(
        file_path=f'{d}/atlas_retrieval_2.json'
    )


def check_file_exists(file_path: str) -> bool:
    """
    Check if the generated file exists.
    Returns True if file exists, False otherwise.
    """
    return os.path.exists(file_path)