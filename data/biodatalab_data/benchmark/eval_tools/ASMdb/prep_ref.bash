#!/bin/bash

# Define input and output file paths
# Input files are paired-end FASTQ reads
R1_IN="data/biodatalab_data/benchmark/results/ASMdb/raw_data/SRR7207779_1.fastq"
R2_IN="data/biodatalab_data/benchmark/results/ASMdb/raw_data/SRR7207779_2.fastq"

# Define the output directory
OUT_DIR="fastp_results"

# Create the output directory if it doesn't exist
mkdir -p "$OUT_DIR"

# Define output file names for trimmed reads
R1_OUT="${OUT_DIR}/SRR7207779_1.trimmed.fastq"
R2_OUT="${OUT_DIR}/SRR7207779_2.trimmed.fastq"

# Define report file names
HTML_REPORT="${OUT_DIR}/SRR7207779.html"
JSON_REPORT="${OUT_DIR}/SRR7207779.json"

# Print the parameters being used
echo "Running Fastp with the following parameters:"
echo "Input Read 1: ${R1_IN}"
echo "Input Read 2: ${R2_IN}"
echo "Output Read 1: ${R1_OUT}"
echo "Output Read 2: ${R2_OUT}"
echo "Sliding window size (-W): 4"
echo "Mean quality requirement (-M): 20"
echo "Quality threshold for a base (-q): 15"
echo "Percentage of unqualified bases allowed (-u): 40"
echo "Max N bases per read (-n): 5"
echo "Low complexity filter threshold (-Y): 0 (disabling the filter)"
echo "HTML Report: ${HTML_REPORT}"
echo "JSON Report: ${JSON_REPORT}"
echo "-------------------------------------------"

# Execute the Fastp command
fastp \
    --in1 "$R1_IN" \
    --in2 "$R2_IN" \
    --out1 "$R1_OUT" \
    --out2 "$R2_OUT" \
    --window_size 4 \
    --qualified_quality_phred 15 \
    --unqualified_percent_limit 40 \
    --n_base_limit 5 \
    --average_qual 20 \
    --low_complexity_filter \
    --complexity_threshold 0 \
    --html "$HTML_REPORT" \
    --json "$JSON_REPORT"

echo "Fastp trimming completed."
echo "Trimmed files are located in the '${OUT_DIR}' directory."