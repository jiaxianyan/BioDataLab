def compodynamics_integration(d):
    try:
        pred_fna_path = f"{d}/compodynamics_integration.fna"
        gold_fna_path = "benchmark/gold_results/compodynamics_integration.fasta"

        def parse_fasta(path):
            sequences = {}
            current_id = None
            current_seq = []

            with open(path, "r") as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    if line.startswith(">"):
                        if current_id is not None:
                            sequences[current_id] = "".join(current_seq).upper()
                        current_id = line[1:].strip()
                        current_seq = []
                    else:
                        current_seq.append(line)

                # save last sequence
                if current_id is not None:
                    sequences[current_id] = "".join(current_seq).upper()

            return sequences

        pred_seqs = parse_fasta(pred_fna_path)
        gold_seqs = parse_fasta(gold_fna_path)

        # 1. same number of sequences
        if len(pred_seqs) != len(gold_seqs):
            return False

        # 2. same sequence IDs
        if set(pred_seqs.keys()) != set(gold_seqs.keys()):
            return False

        # 3. same sequence content
        for seq_id in pred_seqs:
            if pred_seqs[seq_id] != gold_seqs[seq_id]:
                return False

        return True

    except Exception as e:
        print(f"Error processing FNA files: {e}")
        return False
