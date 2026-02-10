def fusionneoantigen_annotate(d):
    """
    Evaluate whether the predicted fusion neoantigen peptide(s)
    match the ground truth sequences.

    Comparison logic:
    - Read FASTA from prediction directory and gold directory
    - Extract sequences only (ignore headers)
    - Order-independent comparison
    """
    try:
        pred_fasta = f"{d}/fusionneoantigen_annotate.fasta"
        gold_fasta = "benchmark/gold_results/fusionneoantigen_annotate.fasta"

        def read_fasta_sequences(path):
            sequences = []
            current_seq = []

            with open(path, "r") as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    if line.startswith(">"):
                        if current_seq:
                            sequences.append("".join(current_seq))
                            current_seq = []
                    else:
                        current_seq.append(line)
                if current_seq:
                    sequences.append("".join(current_seq))

            return sequences

        pred_seqs = read_fasta_sequences(pred_fasta)
        gold_seqs = read_fasta_sequences(gold_fasta)

        # Order-independent comparison
        return sorted(pred_seqs) == sorted(gold_seqs)

    except Exception as e:
        print(f"Error evaluating fusionneoantigen_annotate: {e}")
        return False
