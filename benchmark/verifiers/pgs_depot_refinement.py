import os

def pgs_depot_refinement(d):
    """
    Evaluates the generated VCF file against the ground truth.
    Focuses on CHROM, POS, ID, REF, and ALT columns.
    """
    def get_vcf_variants(filepath):
        variants = set()
        try:
            with open(filepath, 'r') as f:
                for line in f:
                    # Skip metadata headers
                    if line.startswith('##'):
                        continue
                    # Skip the column header line
                    if line.startswith('#CHROM'):
                        continue
                    
                    parts = line.strip().split('\t')
                    if len(parts) >= 5:
                        # Extract: CHROM, POS, ID, REF, ALT
                        # We use a tuple so the data is hashable for the set
                        variant_key = (parts[0], parts[1], parts[2], parts[3], parts[4])
                        variants.add(variant_key)
        except Exception as e:
            print(f"Error reading {filepath}: {e}")
            return None
        return variants

    try:
        # Define paths for predicted and gold standard results
        pred_path = os.path.join(d, 'pgs_depot_refinement.vcf')
        gold_path = 'benchmark/gold_results/pgs_depot_refinement.vcf'

        if not os.path.exists(pred_path):
            print(f"Prediction file not found: {pred_path}")
            return False

        pred_variants = get_vcf_variants(pred_path)
        gold_variants = get_vcf_variants(gold_path)

        if pred_variants is None or gold_variants is None:
            return False

        # Compare the sets of variants
        return pred_variants == gold_variants

    except Exception as e:
        print(f"Evaluation error: {e}")
        return False