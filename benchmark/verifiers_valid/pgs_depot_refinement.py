import os

def pgs_depot_refinement(d):
    """
    Checks if the generated VCF file exists.
    """
    try:
        # Define path for predicted results
        pred_path = os.path.join(d, 'pgs_depot_refinement.vcf')
        
        # Check if file exists
        return os.path.exists(pred_path)
    
    except Exception as e:
        print(f"Evaluation error: {e}")
        return False