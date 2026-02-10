def cancerscem_annotate_2(d):
    try:
        pred_file = f"{d}/cancerscem_annotate_2.txt"
        gold_file = "benchmark/gold_results/cancerscem_annotate_2.txt"

        def load_gene_pairs(path):
            pairs = set()
            with open(path, "r") as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    # 支持空格或 tab 分隔
                    parts = line.split()
                    if len(parts) != 2:
                        continue
                    g1, g2 = parts
                    # 忽略大小写 + 忽略基因顺序
                    pair = tuple(sorted([g1.lower(), g2.lower()]))
                    pairs.add(pair)
            return pairs

        pred_pairs = load_gene_pairs(pred_file)
        gold_pairs = load_gene_pairs(gold_file)

        return pred_pairs == gold_pairs

    except Exception as e:
        print(f"Error processing cancerscem_annotate_2 results: {e}")
        return False
