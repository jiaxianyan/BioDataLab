from rdkit import Chem
import os

mol_dir = '/root/biodatalab/benchmark/dataset/FusionNeoAntigen/dock/mols'
mol_files = os.listdir(mol_dir)

for mol_f in mol_files:
    if mol_f.endswith('pdbqt'):
        continue
    file_path = os.path.join(mol_dir, mol_f)
    mol = Chem.MolFromMol2File(file_path)
    smi = Chem.MolToSmiles(mol)
    
    print(mol_f,smi)