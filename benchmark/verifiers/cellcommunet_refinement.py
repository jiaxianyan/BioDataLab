import os
import anndata
import pandas as pd

def cellcommunet_refinement(d):
    """
    Evaluates the cellcommunet_refinement task.
    
    Args:
        d (str): The directory path containing the generated result file.
                 Example: 'output/result'
    
    Returns:
        bool: True if the generated result matches the ground truth, False otherwise.
    """
    # 1. Define paths
    # Prediction path: <d>/cellcommunet_refinement.h5ad
    pred_path = os.path.join(d, 'cellcommunet_refinement.h5ad')
    
    # Ground Truth path: benchmark/gold_results/cellcommunet_refinement.h5ad
    # (Assuming you follow the same directory structure as your example)
    gold_path = 'benchmark/gold_results/cellcommunet_refinement.h5ad'

    try:
        # 2. Check if prediction file exists
        if not os.path.exists(pred_path):
            print(f"Error: Prediction file not found at {pred_path}")
            return False

        # 3. Load the AnnData objects
        # using 'r' (read-only) mode is safer and faster
        pred_adata = anndata.read_h5ad(pred_path, backed='r')
        gold_adata = anndata.read_h5ad(gold_path, backed='r')

        # 4. Evaluation Logic
        
        # Check 1: Shape comparison (Fastest check)
        # Are the number of cells (n_obs) and genes (n_vars) identical?
        if pred_adata.shape != gold_adata.shape:
            print(f"Shape Mismatch: Pred {pred_adata.shape} vs Gold {gold_adata.shape}")
            return False

        # Check 2: Index comparison (Crucial check)
        # Did the model filter the exact same cells?
        # We sort them to ensure order doesn't affect equality (though usually order is preserved)
        pred_cells = sorted(pred_adata.obs_names.astype(str))
        gold_cells = sorted(gold_adata.obs_names.astype(str))
        
        if pred_cells != gold_cells:
            print("Mismatch in filtered cell IDs (obs_names).")
            return False

        # Check 3 (Optional but strict): Variable/Gene names comparison
        pred_genes = sorted(pred_adata.var_names.astype(str))
        gold_genes = sorted(gold_adata.var_names.astype(str))
        
        if pred_genes != gold_genes:
            print("Mismatch in gene IDs (var_names).")
            return False
            
        # If all checks pass
        return True

    except Exception as e:
        print(f"Error evaluating cellcommunet_refinement: {e}")
        return False