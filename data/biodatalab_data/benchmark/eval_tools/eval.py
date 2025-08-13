#!/usr/bin/env python3
import argparse
import sys

def parse_fasta(file_path):
    """
    Parses a FASTA file and returns a dictionary.

    Args:
        file_path (str): The path to the FASTA file.

    Returns:
        dict: A dictionary where keys are sequence headers (without '>')
              and values are the corresponding sequences.
    """
    sequences = {}
    try:
        with open(file_path, 'r') as f:
            header = None
            current_sequence = []
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # If we have a previous sequence, save it
                    if header:
                        sequences[header] = ''.join(current_sequence)
                    
                    # Start a new sequence
                    header = line[1:] # Store header without the '>'
                    current_sequence = []
                elif header:
                    # Append sequence line
                    current_sequence.append(line)
            
            # Don't forget to save the very last sequence in the file
            if header:
                sequences[header] = ''.join(current_sequence)
                
    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.", file=sys.stderr)
        sys.exit(1)
        
    return sequences

def compare_sequences(dict1, dict2, file1_name, file2_name):
    """
    Compares two sequence dictionaries and prints a detailed report.
    """
    # Python's dictionary comparison is a powerful shortcut!
    # It checks for same keys and same values for each key.
    if dict1 == dict2:
        print("✅ The two FASTA files contain the exact same set of sequences and headers.")
        print(f"   - Found {len(dict1)} sequences in both files.")
        return True

    print("❌ The two FASTA files are different. Here's a summary:\n")
    
    # Get the set of headers from each dictionary
    headers1 = set(dict1.keys())
    headers2 = set(dict2.keys())

    # Find headers unique to each file
    unique_to_1 = headers1 - headers2
    unique_to_2 = headers2 - headers1
    
    # Find headers that are common to both files
    common_headers = headers1 & headers2

    # Check for differences in content for common headers
    differing_content = []
    for header in common_headers:
        if dict1[header] != dict2[header]:
            differing_content.append(header)
    
    # --- Generate Report ---
    is_different = False

    # Report on total sequence counts
    if len(dict1) != len(dict2):
        is_different = True
        print(f"Mismatched sequence counts:")
        print(f"  - '{file1_name}' has {len(dict1)} sequences.")
        print(f"  - '{file2_name}' has {len(dict2)} sequences.\n")

    # Report on headers unique to file 1
    if unique_to_1:
        is_different = True
        print(f"Headers found only in '{file1_name}' ({len(unique_to_1)}):")
        # Print first 5 for brevity
        for header in list(unique_to_1)[:5]:
            print(f"  - {header}")
        if len(unique_to_1) > 5:
            print(f"  - ... and {len(unique_to_1) - 5} more.")
        print()

    # Report on headers unique to file 2
    if unique_to_2:
        is_different = True
        print(f"Headers found only in '{file2_name}' ({len(unique_to_2)}):")
        # Print first 5 for brevity
        for header in list(unique_to_2)[:5]:
            print(f"  - {header}")
        if len(unique_to_2) > 5:
            print(f"  - ... and {len(unique_to_2) - 5} more.")
        print()

    # Report on sequences with same header but different content
    if differing_content:
        is_different = True
        print(f"Headers with different sequences ({len(differing_content)}):")
        # Print first 5 for brevity
        for header in differing_content[:5]:
            print(f"  - {header}")
        if len(differing_content) > 5:
            print(f"  - ... and {len(differing_content) - 5} more.")
        print()

    return not is_different

def eval_1(res_data, ref_data):
    """Main function to run the script."""
    parser = argparse.ArgumentParser(
        description="Compare two FASTA files to see if they contain the same sequences.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    args = parser.parse_args()
    
    print(f"Comparing '{res_data}' and '{ref_data}'...\n")

    # Step 1: Parse both files into dictionaries
    seqs1 = parse_fasta(res_data)
    seqs2 = parse_fasta(ref_data)
    
    # Step 2: Compare the dictionaries and report the results
    return compare_sequences(seqs1, seqs2, res_data, ref_data)

if __name__ == '__main__':
    eval_1('data/biodatalab_data/benchmark/results/1/ref.fasta',
           'data/biodatalab_data/benchmark/results/1/ref.fasta')