
def prepare_1(input_fasta, ref_fasta):
    from Bio import SeqIO

    # Define length constraints
    min_len = 16
    max_len = 2700

    # Create a list to hold the filtered sequences
    filtered_sequences = []

    print(f"Starting to filter {input_fasta}...")
    print(f"Keeping sequences with length between {min_len} and {max_len} residues.")

    # Iterate through the sequences and filter them
    for record in SeqIO.parse(input_fasta, "fasta"):
        if min_len <= len(record.seq) <= max_len:
            filtered_sequences.append(record)

    # Write the filtered sequences to the output file
    SeqIO.write(filtered_sequences, ref_fasta, "fasta")

    print(f"\nFiltering complete.")
    print(f"A total of {len(filtered_sequences)} sequences were saved to {ref_fasta}.")

if __name__ == '__main__':
    prepare_1('./data/biodatalab_data/data_lake/uniprot_sprot.fasta', 'data/biodatalab_data/benchmark/results/1/ref.fasta')