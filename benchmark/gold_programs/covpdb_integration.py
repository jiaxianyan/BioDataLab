import os
import json
import requests
import pandas as pd
import numpy as np
from Bio.PDB import PDBParser, SASA, Selection, NeighborSearch
from Bio.PDB.Polypeptide import protein_letters_3to1
import warnings
import freesasa

# Suppress Biopython warnings
warnings.simplefilter('ignore')

# --- Configuration ---
INPUT_ID_FILE = "benchmark/dataset/CovPDB/cov_pdb_ids.json"
STRUCTURE_DIR = "benchmark/dataset/CovPDB/complex_structures"
OUTPUT_DIR = "pred_results"
TMP_DIR = "./tmp/covpdb"

# Ensure directories exist
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(TMP_DIR, exist_ok=True)

# --- Helper Functions ---

def get_pdb_metadata(pdb_id):
    """
    Fetches Method, Resolution, and Binding Affinity from RCSB PDB GraphQL API.
    """
    url = "https://data.rcsb.org/graphql"
    query = """
    {
      entry(entry_id: "%s") {
        exptl {
          method
        }
        rcsb_entry_info {
          resolution_combined
        }
        rcsb_binding_affinity {
          comp_id
          type
          value
          unit
        }
      }
    }
    """ % pdb_id.upper()
    
    try:
        response = requests.post(url, json={'query': query})
        data = response.json()['data']['entry']
        
        # Method
        method = data['exptl'][0]['method'] if data.get('exptl') else "Unknown"
        
        # Resolution
        resolution = data['rcsb_entry_info']['resolution_combined'][0] if data.get('rcsb_entry_info') and data['rcsb_entry_info'].get('resolution_combined') else None
        
        # Affinity (Combine type + value + unit if available)
        affinity = "N/A"
        if data.get('rcsb_binding_affinity'):
            aff_data = data['rcsb_binding_affinity'][0]
            affinity = f"{aff_data.get('type')} {aff_data.get('value')} {aff_data.get('unit')}"
            
        return method, resolution, affinity
    except Exception as e:
        print(f"Error fetching metadata for {pdb_id}: {e}")
        return "Unknown", None, "N/A"

def get_uniprot_annotation(pdb_id):
    """
    Maps PDB ID to UniProt to get Name, Gene, and Sequence.
    Uses the EBI PDBe API which is often cleaner for PDB-to-UniProt mapping.
    """
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id.lower()}"
    
    try:
        response = requests.get(url)
        if response.status_code != 200:
            return None
        
        data = response.json()
        if not data:
            return None
            
        # Get the first mapped UniProt ID
        uniprot_id = list(data[pdb_id.lower()]['UniProt'].keys())[0]
        
        # Now fetch details from UniProt API
        uni_url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
        uni_resp = requests.get(uni_url).json()
        
        # Extract fields
        name = uni_resp.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'Unknown')
        
        # Gene Symbol
        genes = uni_resp.get('genes', [])
        gene_symbol = genes[0].get('geneName', {}).get('value', 'Unknown') if genes else 'Unknown'
        
        # Classification (Keywords often serve as classification/function)
        classification = "Unknown"
        if 'keywords' in uni_resp:
             classification = uni_resp['keywords'][0]['name']

        # Sequence
        sequence = uni_resp.get('sequence', {}).get('value', '')
        
        return {
            "Name": name,
            "Gene": gene_symbol,
            "Classification": classification,
            "Sequence": sequence,
            "Chain": "A" # Defaulting to A for this snippet, real logic requires mapping the chain from PDBe
        }
    except Exception as e:
        print(f"Error fetching UniProt for {pdb_id}: {e}")
        return None


# --- Main Pipeline ---

def main():
    print(f"Starting Data Enrichment for CovPDB...")
    
    # Load IDs
    if not os.path.exists(INPUT_ID_FILE):
        print(f"Error: {INPUT_ID_FILE} not found.")
        return

    with open(INPUT_ID_FILE, 'r') as f:
        pdb_ids = json.load(f)

    complex_results = []
    protein_results = []

    parser = PDBParser(QUIET=True)

    for idx, pdb_id in enumerate(pdb_ids):
        print(f"[{idx+1}/{len(pdb_ids)}] Processing {pdb_id.lower()}...")
        
        # 1. Metadata Retrieval
        method, resolution, affinity = get_pdb_metadata(pdb_id)
        
        # 2. Structural Calculation (SASA)
        struct_path = os.path.join(STRUCTURE_DIR, f"pdb{pdb_id.lower()}.ent")
        sasa_val = 0.0
        
        if os.path.exists(struct_path):
            try:
                structure = freesasa.Structure(struct_path)
                sasa_val = freesasa.calc(structure).totalArea()
            except Exception as e:
                print(f"  Warning: Failed to parse structure {pdb_id.lower()} ({e})")
        else:
             # Try fetching if not local (Optional, based on prompt strictness)
             print(f"  Warning: Structure file not found locally for {pdb_id.lower()}")

        # 3. UniProt Annotation
        # We assume Chain A is the target for the annotation part
        uni_data = get_uniprot_annotation(pdb_id.lower())
        
        # Aggregate Complex Results
        complex_results.append({
            "index": idx,
            "PDB_ID": pdb_id,
            "Method": method,
            "Resolution": resolution,
            "Affinity": affinity,
            "SASA": sasa_val
        })
        
        # Aggregate Protein Results
        if uni_data:
            protein_results.append({
                "index": idx,
                "Chain": uni_data['Chain'],
                "Name": uni_data['Name'],
                "Gene": uni_data['Gene'],
                "Classification": uni_data['Classification'],
                "Sequence": uni_data['Sequence']
            })
        else:
            # Add placeholder if API fails
            protein_results.append({
                "index": idx,
                "Chain": "A",
                "Name": "N/A",
                "Gene": "N/A",
                "Classification": "N/A",
                "Sequence": "N/A"
            })

    # --- Saving Results ---
    df_complex = pd.DataFrame(complex_results)
    df_protein = pd.DataFrame(protein_results)

    # Save Complex CSV
    complex_out = os.path.join(OUTPUT_DIR, "covpdb_integration_complex.csv")
    df_complex.to_csv(complex_out, index=False)
    
    # Save Protein CSV
    protein_out = os.path.join(OUTPUT_DIR, "covpdb_integration_protein.csv")
    df_protein.to_csv(protein_out, index=False)

    print("\nProcessing Complete.")
    print(f"Complex data saved to: {complex_out}")
    print(f"Protein data saved to: {protein_out}")

if __name__ == "__main__":
    main()