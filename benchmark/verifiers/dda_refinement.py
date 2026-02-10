import gzip

def dda_refinement(d):
    try:
        pred_path = f"{d}/dda_refinement.fastq.gz"
        gold_path = "benchmark/gold_results/dda_refinement.fastq.gz"

        def load_sequences(fastq_gz):
            seqs = []
            with gzip.open(fastq_gz, "rt") as f:
                while True:
                    header = f.readline()
                    if not header:
                        break
                    seq = f.readline().strip()
                    f.readline()  # +
                    f.readline()  # quality
                    seqs.append(seq)
            return seqs

        pred_seqs = load_sequences(pred_path)
        gold_seqs = load_sequences(gold_path)

        # 顺序不敏感比较
        return sorted(pred_seqs) == sorted(gold_seqs)

    except Exception as e:
        print(f"Error processing FASTQ files: {e}")
        return False
