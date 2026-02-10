import numpy as np

def get_mol2_center(file_path):
    coords = []
    with open(file_path, 'r') as f:
        reading_atoms = False
        for line in f:
            if "@<TRIPOS>ATOM" in line:
                reading_atoms = True
                continue
            if "@<TRIPOS>BOND" in line:
                reading_atoms = False
                break
            
            if reading_atoms:
                parts = line.split()
                if len(parts) >= 5:
                    # MOL2 格式中坐标通常在第 3, 4, 5 列 (index 2, 3, 4)
                    x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                    coords.append([x, y, z])
    
    coords = np.array(coords)
    center = np.mean(coords, axis=0)
    # 计算范围（用于参考 size）
    size = np.max(coords, axis=0) - np.min(coords, axis=0)
    
    return center, size

ligand_path = "/root/biodatalab/benchmark/dataset/FusionNeoAntigen/dock/1a1e_ligand.mol2"
center, size = get_mol2_center(ligand_path)

print(f"--- Vina Configuration Parameters ---")
print(f"center_x = {center[0]:.3f}")
print(f"center_y = {center[1]:.3f}")
print(f"center_z = {center[2]:.3f}")
# 通常在原始分子尺寸基础上各边加 10-15埃 的 buffer
print(f"size_x = {size[0] + 12:.3f}")
print(f"size_y = {size[1] + 12:.3f}")
print(f"size_z = {size[2] + 12:.3f}")