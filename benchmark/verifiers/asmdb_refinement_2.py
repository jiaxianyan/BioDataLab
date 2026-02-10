import os
import pysam

def asmdb_refinement_2(d):
    """
    Evaluate if the generated BAM file is biologically equivalent to the gold standard.
    
    Args:
        d (str): The directory containing the generated 'asmdb_refinement_2.bam'.
    """
    pred_bam_path = os.path.join(d, 'asmdb_refinement_2.bam')
    gold_bam_path = 'benchmark/gold_results/asmdb_refinement_2.bam'

    # 1. Check if the prediction file exists
    if not os.path.exists(pred_bam_path):
        print(f"Error: Output file not found at {pred_bam_path}")
        return False

    try:
        # Open both BAM files using pysam
        # 'rb' means read binary (BAM format)
        with pysam.AlignmentFile(pred_bam_path, "rb") as pred_bam, \
             pysam.AlignmentFile(gold_bam_path, "rb") as gold_bam:

            # --- Check 1: Verify Reference Genome Consistency ---
            # Checks if the chromosome names and lengths are identical
            if pred_bam.header.references != gold_bam.header.references:
                print("Mismatch: Reference genome headers (chromosomes) do not match.")
                return False
            
            if pred_bam.header.lengths != gold_bam.header.lengths:
                print("Mismatch: Reference genome lengths do not match.")
                return False

            # --- Check 2: Verify Alignment Statistics ---
            # .mapped: Number of mapped reads
            # .unmapped: Number of unmapped reads
            # This is robust against timestamp/header ID differences
            
            if pred_bam.mapped != gold_bam.mapped:
                print(f"Mismatch: Mapped reads count. Pred: {pred_bam.mapped}, Gold: {gold_bam.mapped}")
                return False
            
            if pred_bam.unmapped != gold_bam.unmapped:
                print(f"Mismatch: Unmapped reads count. Pred: {pred_bam.unmapped}, Gold: {gold_bam.unmapped}")
                return False

            # (Optional Strict Check) Verify coordinate sorting
            # Check if the file is marked as sorted in the header
            if 'HD' in pred_bam.header and 'SO' in pred_bam.header['HD']:
                 if pred_bam.header['HD']['SO'] != 'coordinate':
                     print("Mismatch: BAM file is not sorted by coordinate.")
                     return False

            return True

    except ValueError as e:
        print(f"Error: File is not a valid BAM or is corrupt. {e}")
        return False
    except Exception as e:
        print(f"An unexpected error occurred during evaluation: {e}")
        return False