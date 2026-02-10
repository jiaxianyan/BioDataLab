def asmdb_refinement_1(d):
    try:
        pred_path = f"{d}/asmdb_refinement_1.bed"
        gold_path = "benchmark/gold_results/asmdb_refinement_1.bed"

        def load_bed_sites(path):
            sites = set()
            with open(path, "r") as f:
                for line in f:
                    if line.strip() == "" or line.startswith("#"):
                        continue
                    fields = line.strip().split("\t")
                    # BED core fields we care about
                    chrom = fields[0]
                    start = int(fields[1])
                    end = int(fields[2])
                    strand = fields[5]
                    sites.add((chrom, start, end, strand))
            return sites

        pred_sites = load_bed_sites(pred_path)
        gold_sites = load_bed_sites(gold_path)

        return pred_sites == gold_sites

    except Exception as e:
        print(f"Error evaluating asmdb_refinement_1: {e}")
        return False
