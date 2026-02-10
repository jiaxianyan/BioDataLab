import sys

print("Locus_Tag\tSymbol\tOld_Locus_Tag\tLocation\tProtein_ID\tProduct")

num = 0
with open("benchmark/dataset/CyanoOmicsDB/ncbi_dataset/data/demo_genomic.gff", "r") as f:
    for line in f:
        if line.startswith("#") or "\tCDS\t" not in line:
            continue
        num += 1
        cols = line.split("\t")
        # 提取基因组位置 (第1, 4, 5, 7列: 染色体, 起始, 终止, 正负链)
        location = f"{cols[0]}:{cols[3]}-{cols[4]}({cols[6]})"
        
        # 提取第9列的属性
        attrs = {item.split("=")[0]: item.split("=")[1] for item in cols[8].strip().split(";") if "=" in item}
        
        print(f"{attrs.get('locus_tag', '-')}\t"
              f"{attrs.get('gene', '-')}\t"
              f"{attrs.get('old_locus_tag', '-')}\t"
              f"{location}\t"
              f"{attrs.get('protein_id', '-')}\t"
              f"{attrs.get('product', '-')}")
        
print(num)