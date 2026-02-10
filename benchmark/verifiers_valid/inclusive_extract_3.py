import os


def inclusive_extract_3(d) -> bool:
    return check_file_exists(
        out_json_path = f"{d}/inclusive_extract_3.json",
    )

def check_file_exists(out_json_path: str) -> bool:
    """
    Check if the generated file exists.
    
    Args:
        out_json_path: Path to the generated JSON file
    
    Returns:
        True if file exists, False otherwise
    """
    return os.path.exists(out_json_path)