import os

def asmdb_refinement_2(d):
    """
    Check if the generated BAM file exists.
    
    Args:
        d (str): The directory containing the generated 'asmdb_refinement_2.bam'.
    
    Returns:
        bool: True if the file exists, False otherwise.
    """
    pred_bam_path = os.path.join(d, 'asmdb_refinement_2.bam')
    
    return os.path.exists(pred_bam_path)