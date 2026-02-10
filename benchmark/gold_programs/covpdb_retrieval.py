import os
import json
import requests
import warnings
from Bio.PDB import PDBList, MMCIFParser, NeighborSearch, Selection
from Bio.PDB.PDBExceptions import PDBConstructionWarning

# Suppress Biopython warnings for cleaner output
warnings.simplefilter('ignore', PDBConstructionWarning)

# --- Configuration ---
START_DATE = "2018-10-01T00:00:00Z"
END_DATE = "2018-10-07T00:00:00Z"
RES_THRESHOLD = 2.5
COVALENT_DIST_THRESHOLD = 2.0
ARTIFACT_FILE_PATH = "benchmark/dataset/CovPDB/ligand_list"
OUTPUT_DIR = "pred_results"
OUTPUT_FILE = os.path.join(OUTPUT_DIR, "cov_pdb_retrieval.json")
TEMP_PDB_DIR = "temp_pdb_structures"

def get_pdbs_from_rcsb():
    """
    Queries RCSB PDB Search API for entries within the date range
    and resolution threshold.
    """
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_accession_info.initial_release_date",
                        "operator": "range",
                        "value": {
                            "from": START_DATE,
                            "to": END_DATE
                        }                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.resolution_combined",
                        "operator": "less_or_equal",
                        "value": RES_THRESHOLD
                    }
                }
            ]
        },
        "request_options": {
            "return_all_hits": True
        },
        "return_type": "entry"
    }

    try:
        response = requests.post(url, json=query)
        response.raise_for_status()
        data = response.json()
        pdb_ids = [entry['identifier'] for entry in data.get('result_set', [])]
        print(f"INFO: Found {len(pdb_ids)} candidate PDBs from RCSB Search.")
        return pdb_ids
    except Exception as e:
        print(f"ERROR: Failed to query RCSB API. {e}")
        return []

def load_artifact_list(filepath):
    """
    Loads the BioLiP artifact list. Assumes a text file with one ligand code per line
    or space-separated codes.
    """
    artifacts = set()
    if not os.path.exists(filepath):
        print(f"WARNING: Artifact file not found at {filepath}. No ligands will be filtered as artifacts.")
        return artifacts
    
    try:
        with open(filepath, 'r') as f:
            content = f.read().strip().split('\n')
            # Split by whitespace or newlines to handle various formats
            for c in content:
                artifacts.add(c.split()[0].upper())
    except Exception as e:
        print(f"ERROR: Could not read artifact file. {e}")
    
    return artifacts

def is_covalent_complex(pdb_id, artifact_list):
    """
    Downloads PDB, parses it, and checks for covalent bonding (< 2.0 A)
    between protein and non-artifact ligands.
    """
    pdbl = PDBList(verbose=False)
    
    try:
        # Download mmCIF file
        cif_file = pdbl.retrieve_pdb_file(pdb_id, pdir=TEMP_PDB_DIR, file_format='pdb')
        
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure(pdb_id, cif_file)
        
        # Get all atoms and categorize
        atom_list = [atom for atom in structure.get_atoms()]
        ns = NeighborSearch(atom_list)
        
        # Iterate over all residues to find potential ligands
        for residue in structure.get_residues():
            # Filter 1: Check if it is a HETATM (ligand/water)
            # 'H_ ' is HETATM, 'W' is water in Biopython id[0]
            het_flag = residue.id[0].strip()
            res_name = residue.resname.upper()

            # Skip standard ATOM records (blank het_flag) and Water
            if het_flag == '' or 'W' in het_flag or res_name == 'HOH':
                continue
            
            # Filter 2: Skip Artifacts
            if res_name in artifact_list:
                continue

            # Check distance to Protein
            # We look for any protein atom within COVALENT_DIST_THRESHOLD of any ligand atom
            ligand_atoms = list(residue.get_atoms())
            
            for lig_atom in ligand_atoms:
                # Search neighbors < 2.0 A
                neighbors = ns.search(lig_atom.get_coord(), COVALENT_DIST_THRESHOLD)
                
                for neighbor in neighbors:
                    # Check if neighbor is part of the Protein (standard residues)
                    neighbor_res = neighbor.get_parent()
                    neighbor_flag = neighbor_res.id[0].strip()
                    
                    # If neighbor is a standard amino acid (flag is empty) and not the ligand itself
                    if neighbor_flag == '' and neighbor_res != residue:
                        # Found a covalent bond!
                        # Optional: Print detail
                        # print(f"  Hit: {pdb_id} {res_name} bound to {neighbor_res.resname}")
                        return True
                        
    except Exception as e:
        print(f"ERROR: Processing {pdb_id} failed: {e}")
        return False
        
    return False

def main():
    # 1. Setup directories
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(TEMP_PDB_DIR, exist_ok=True)
    
    # 2. Get Candidates
    candidate_pdbs = get_pdbs_from_rcsb()
    if not candidate_pdbs:
        return

    # 3. Load Artifacts
    artifact_list = load_artifact_list(ARTIFACT_FILE_PATH)
    print(f"INFO: Loaded {len(artifact_list)} artifact ligands.")

    # 4. Filter for Covalent Bonding
    covalent_pdbs = []
    print("INFO: Starting geometric validation (this may take a moment)...")
    
    for i, pdb_id in enumerate(candidate_pdbs):
        if is_covalent_complex(pdb_id, artifact_list):
            covalent_pdbs.append(pdb_id)
        
        # Simple progress logger
        if (i + 1) % 10 == 0:
            print(f"Processed {i + 1}/{len(candidate_pdbs)}...")

    # 5. Save Results
    print(f"INFO: Identified {len(covalent_pdbs)} covalent complexes.")
    with open(OUTPUT_FILE, 'w') as f:
        json.dump(covalent_pdbs, f, indent=4)
    
    print(f"SUCCESS: Results saved to {OUTPUT_FILE}")

    # Optional: Cleanup temp directory
    # import shutil
    # shutil.rmtree(TEMP_PDB_DIR)

if __name__ == "__main__":
    main()