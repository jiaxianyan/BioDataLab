import mdtraj as md
import numpy as np
import json
import os
import subprocess
import sys

# ==========================================
# Configuration & Path Setup
# ==========================================
TRAJ_PATH = "benchmark/dataset/SCoV2-MD/474_trj_245.xtc"
TOPO_PATH = "benchmark/dataset/SCoV2-MD/473_dyn_245.pdb"
OUTPUT_DIR = "benchmark/gold_results"

# Ensure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

def calculate_rmsf(traj, output_path):
    """
    Calculates RMSF for Alpha-Carbon (CA) atoms after superposing 
    the trajectory to the first frame to remove rigid body motion.
    """
    print("Starting RMSF Calculation...")
    
    # 1. Select Alpha-Carbon (CA) atoms
    ca_indices = traj.topology.select('name CA')
    
    # 2. Slice trajectory to only include CA atoms
    traj_ca = traj.atom_slice(ca_indices)
    
    # 3. Superpose trajectory to the first frame (reference) 
    # using the CA atoms to align.
    traj_ca.superpose(traj_ca, 0)
    
    # 4. Calculate RMSF (relative to the first frame after alignment)
    # returns an array of shape (n_atoms,)
    rmsf_values = md.rmsf(traj_ca, traj_ca, 0)
    
    # 5. Format results: 0-based residue index -> RMSF (nm)
    # Since we sliced CA atoms, the list index corresponds to the residue count.
    results = {
        i: float(value) 
        for i, value in enumerate(rmsf_values)
    }
    
    # 6. Save to JSON
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=4)
        
    print(f"RMSF results saved to: {output_path}")


def main():
    # Load Trajectory
    print(f"Loading trajectory: {TRAJ_PATH}")
    try:
        traj = md.load(TRAJ_PATH, top=TOPO_PATH)
    except IOError as e:
        print(f"Error loading files: {e}")
        sys.exit(1)
        
    # Task 1: RMSF
    rmsf_out = os.path.join(OUTPUT_DIR, "rmsf_ca.json")
    calculate_rmsf(traj, rmsf_out)

    
    print("\nBenchmark generation complete.")

if __name__ == "__main__":
    main()