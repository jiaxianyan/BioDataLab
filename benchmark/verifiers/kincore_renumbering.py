from Bio.PDB import PDBParser
import numpy as np
from typing import Dict, Tuple


def kincore_renumbering(d) -> bool:
    is_match, metrics = compare_pdb_files(
        ref_pdb_path = 'benchmark/gold_results/kincore_renumbering.pdb', 
        out_pdb_path = f"{d}/kincore_renumbering.pdb",
    )
    return is_match
    

def compare_pdb_files(ref_pdb_path: str, out_pdb_path: str) -> Tuple[bool, Dict]:
    """
    Compare two PDB files atom-by-atom using chain id, residue id, residue name,
    insertion code, and atom name as keys, then compare x,y,z coordinates.
    Returns (is_match, metrics) where is_match is True if F1 score == 1.
    """
    try:
        parser = PDBParser()
        ref_structure = parser.get_structure('ref', ref_pdb_path)
        out_structure = parser.get_structure('out', out_pdb_path)
        
        # Create dictionaries to store atoms with unique keys
        def create_atom_dict(structure):
            atom_dict = {}
            for model in structure:
                for chain in model:
                    chain_id = chain.id
                    for residue in chain:
                        res_id = residue.id  # (hetflag, sequence identifier, insertion code)
                        res_name = residue.resname
                        for atom in residue:
                            atom_name = atom.name
                            # Create unique key
                            key = (chain_id, res_id[1], res_name, res_id[2], atom_name)
                            atom_dict[key] = atom.coord
            return atom_dict
        
        ref_atoms = create_atom_dict(ref_structure)
        out_atoms = create_atom_dict(out_structure)
        
        # Calculate matches
        true_positives = 0
        false_positives = 0
        false_negatives = 0
        
        # Check output atoms against reference
        for key, out_coord in out_atoms.items():
            if key in ref_atoms:
                ref_coord = ref_atoms[key]
                # Compare coordinates with tolerance for floating point
                if np.allclose(out_coord, ref_coord, atol=1e-3):
                    true_positives += 1
                else:
                    false_positives += 1
            else:
                false_positives += 1
        
        # Check reference atoms not in output
        for key in ref_atoms:
            if key not in out_atoms:
                false_negatives += 1
        
        # Calculate metrics
        precision = true_positives / (true_positives + false_positives) if (true_positives + false_positives) > 0 else 0
        recall = true_positives / (true_positives + false_negatives) if (true_positives + false_negatives) > 0 else 0
        
        if precision + recall == 0:
            f1_score = 0
        else:
            f1_score = 2 * precision * recall / (precision + recall)
        
        # Check if F1 score equals 1
        is_match = abs(f1_score - 1.0) < 1e-10
        
        metrics = {
            'f1_score': f1_score,
            'precision': precision,
            'recall': recall,
            'true_positives': true_positives,
            'false_positives': false_positives,
            'false_negatives': false_negatives,
            'total_ref_atoms': len(ref_atoms),
            'total_out_atoms': len(out_atoms),
            'is_match': is_match
        }
        
        return is_match, metrics
        
    except Exception as e:
        return False, {'error': str(e)}