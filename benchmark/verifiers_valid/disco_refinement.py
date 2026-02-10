import os

def disco_refinement(d):
    """
    Checks whether the generated RDS file exists.
    
    Args:
        d (str): The directory containing the generated 'disco_refinement.rds'.
                 Example: 'benchmark/output/task_1'
                 
    Returns:
        bool: True if the file exists, False otherwise.
    """
    # Construct path to the generated file
    pred_path = os.path.join(d, 'disco_refinement.rds')
    
    # Check if file exists
    if os.path.exists(pred_path):
        return True
    else:
        return False