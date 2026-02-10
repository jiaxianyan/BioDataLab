from Bio import SeqIO
from typing import Dict, Tuple

def plabdab_retrieval(d) -> bool:
    is_match, metrics = compare_fasta_files(
        ref_fasta_path = 'benchmark/gold_results/antibody_seq_retrieval.fasta', 
        output_fasta_path = f"{d}/antibody_seq_retrieval.fasta", 
        f1_threshold = 0.9,
    )
    return is_match

def compare_fasta_files(
    ref_fasta_path: str, 
    output_fasta_path: str, 
    f1_threshold: float = 0.9
) -> Tuple[bool, Dict]:
    """
    Compare two FASTA files by matching IDs and sequences.
    
    Args:
        ref_fasta_path: Path to reference FASTA file
        output_fasta_path: Path to output FASTA file
        f1_threshold: Threshold for F1 score to consider as match (default: 0.9)
    
    Returns:
        Tuple of (is_match, metrics_dict)
        is_match: True if F1 score >= threshold, False otherwise
        metrics_dict: Dictionary containing detailed comparison metrics
    """
    try:
        # Parse FASTA files into dictionaries
        def parse_fasta_to_dict(fasta_path: str) -> Dict[str, str]:
            """Parse FASTA file into {id: sequence} dictionary."""
            seq_dict = {}
            for record in SeqIO.parse(fasta_path, "fasta"):
                # Use the first word of the ID (before whitespace) as key
                seq_id = record.id
                seq_dict[seq_id] = str(record.seq).upper()
            return seq_dict
        
        # Load sequences
        ref_seqs = parse_fasta_to_dict(ref_fasta_path)
        output_seqs = parse_fasta_to_dict(output_fasta_path)
        
        # Calculate matches
        true_positives = 0
        false_positives = 0
        false_negatives = 0
        
        matched_ids = []
        mismatched_ids = []
        missing_ids = []
        
        # Check output sequences against reference
        for seq_id, output_seq in output_seqs.items():
            if seq_id in ref_seqs:
                if output_seq == ref_seqs[seq_id]:
                    true_positives += 1
                    matched_ids.append(seq_id)
                else:
                    false_positives += 1
                    mismatched_ids.append((seq_id, "sequence_mismatch"))
            else:
                false_positives += 1
                mismatched_ids.append((seq_id, "id_not_in_reference"))
        
        # Check reference sequences not in output
        for seq_id in ref_seqs:
            if seq_id not in output_seqs:
                false_negatives += 1
                missing_ids.append(seq_id)
        
        # Calculate metrics
        precision = true_positives / (true_positives + false_positives) if (true_positives + false_positives) > 0 else 0
        recall = true_positives / (true_positives + false_negatives) if (true_positives + false_negatives) > 0 else 0
        
        if precision + recall == 0:
            f1_score = 0
        else:
            f1_score = 2 * precision * recall / (precision + recall)
        
        # Prepare results
        is_match = f1_score >= f1_threshold
        
        metrics = {
            "f1_score": f1_score,
            "precision": precision,
            "recall": recall,
            "true_positives": true_positives,
            "false_positives": false_positives,
            "false_negatives": false_negatives,
            "total_ref_sequences": len(ref_seqs),
            "total_output_sequences": len(output_seqs),
            "matched_count": len(matched_ids),
            "mismatched_count": len(mismatched_ids),
            "missing_count": len(missing_ids),
            "f1_threshold": f1_threshold,
            "is_match": is_match,
            "matched_ids": matched_ids[:50],  # Limit to first 50 for readability
            "missing_ids": missing_ids[:50],   # Limit to first 50 for readability
            "mismatched_details": mismatched_ids[:50]  # Limit to first 50
        }
        
        return is_match, metrics
        
    except FileNotFoundError as e:
        return False, {"error": f"File not found: {str(e)}"}
    except Exception as e:
        return False, {"error": f"Error processing FASTA files: {str(e)}"}