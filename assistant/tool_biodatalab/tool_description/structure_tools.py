description = [
    {
        "description": "Import FASTQ data via a QIIME 2 manifest, optionally join paired-end reads, "
                       "trim adapters, perform quality-score filtering, and generate demux summary.",
        "name": "q2_ingest_and_qc_api",
        "optional_parameters": [
            {"default": None, "description": "Sample metadata TSV for later visualizations/analyses", "name": "sample_metadata_tsv", "type": "str"},
            {"default": True, "description": "Whether reads are paired-end", "name": "paired_end", "type": "bool"},
            {"default": False, "description": "If paired_end, join pairs using VSEARCH", "name": "join_paired", "type": "bool"},
            {"default": None, "description": "Extra cutadapt parameters dict (e.g. {'p_front_f': 'ACTG'})", "name": "cutadapt_params", "type": "dict"},
            {"default": False, "description": "Enable quality-filter q-score", "name": "qscore_filter", "type": "bool"},
            {"default": None, "description": "Forward truncation length hint (recorded in provenance)", "name": "trunc_len_f", "type": "int"},
            {"default": None, "description": "Reverse truncation length hint (recorded in provenance)", "name": "trunc_len_r", "type": "int"},
            {"default": 33, "description": "Phred offset used in manifest format name (33 or 64)", "name": "phred_offset", "type": "int"},
            {"default": False, "description": "Re-run even if outputs exist", "name": "force", "type": "bool"},
        ],
        "required_parameters": [
            {"default": None, "description": "Path to QIIME 2 manifest CSV (single or paired)", "name": "manifest_csv", "type": "str"},
            {"default": "./q2_ingest_and_qc_out", "description": "Directory to save outputs", "name": "output_dir", "type": "str"},
        ],
    },
    {
        "description": "Denoise reads (Deblur or DADA2), generate feature table and representative sequences, "
                       "and optionally filter samples by minimal reads.",
        "name": "q2_denoise_and_feature_table_api",
        "optional_parameters": [
            {"default": "deblur", "description": "Denoising method: 'deblur' or 'dada2'", "name": "method", "type": "str"},
            {"default": 250, "description": "Trim length for Deblur (ignored for DADA2)", "name": "trim_length", "type": "int"},
            {"default": None, "description": "DADA2 forward trunc length (paired or single)", "name": "trunc_len_f", "type": "int"},
            {"default": None, "description": "DADA2 reverse trunc length (paired)", "name": "trunc_len_r", "type": "int"},
            {"default": None, "description": "Filter out samples with fewer total reads than this", "name": "min_reads_per_sample", "type": "int"},
            {"default": 0, "description": "Threads for denoising (0 lets QIIME decide)", "name": "threads", "type": "int"},
            {"default": False, "description": "Re-run even if outputs exist", "name": "force", "type": "bool"},
        ],
        "required_parameters": [
            {"default": None, "description": "Path to demultiplexed sequences artifact (.qza)", "name": "demux_qza", "type": "str"},
            {"default": "./q2_denoise_out", "description": "Directory to save outputs", "name": "output_dir", "type": "str"},
        ],
    },
    {
        "description": "Assign taxonomy using a pretrained classifier (preferred) or by training from reference, "
                       "then generate taxa barplot visualization.",
        "name": "q2_taxonomy_api",
        "optional_parameters": [
            {"default": None, "description": "Sample metadata TSV for taxa barplot", "name": "metadata_tsv", "type": "str"},
            {"default": None, "description": "Path to pretrained classifier artifact (.qza)", "name": "pretrained_classifier_qza", "type": "str"},
            {"default": None, "description": "Reference sequences artifact (.qza) if training classifier", "name": "ref_seqs_qza", "type": "str"},
            {"default": None, "description": "Reference taxonomy artifact (.qza) if training classifier", "name": "ref_taxonomy_qza", "type": "str"},
            {"default": "sklearn", "description": "Taxonomy assignment method: 'sklearn' or 'vsearch'", "name": "method", "type": "str"},
            {"default": 0.97, "description": "Percent identity threshold for vsearch", "name": "perc_identity", "type": "float"},
            {"default": False, "description": "Re-run even if outputs exist", "name": "force", "type": "bool"},
        ],
        "required_parameters": [
            {"default": None, "description": "Feature table artifact (.qza)", "name": "table_qza", "type": "str"},
            {"default": None, "description": "Representative sequences artifact (.qza)", "name": "rep_seqs_qza", "type": "str"},
            {"default": "./q2_taxonomy_out", "description": "Directory to save outputs", "name": "output_dir", "type": "str"},
        ],
    },
    {
        "description": "Build multiple sequence alignment, mask hypervariable regions, "
                       "construct phylogenetic trees, and compute core diversity metrics (alpha, beta, PCoA, Emperor).",
        "name": "q2_phylogeny_and_core_diversity_api",
        "optional_parameters": [
            {"default": 1000, "description": "Sampling depth for diversity analysis", "name": "sampling_depth", "type": "int"},
            {"default": True, "description": "Use midpoint rooting (if False, use default rooting)", "name": "use_midpoint_root", "type": "bool"},
            {"default": False, "description": "Re-run even if outputs exist", "name": "force", "type": "bool"},
        ],
        "required_parameters": [
            {"default": None, "description": "Feature table artifact (.qza)", "name": "table_qza", "type": "str"},
            {"default": None, "description": "Representative sequences artifact (.qza)", "name": "rep_seqs_qza", "type": "str"},
            {"default": "./q2_phylogeny_out", "description": "Directory to save outputs", "name": "output_dir", "type": "str"},
        ],
    },
    {
        "description": "Run common differential abundance analyses (ANCOM, ALDEx2, Songbird placeholders), "
                       "and export feature table and taxonomy artifacts to BIOM and TSV.",
        "name": "q2_stats_and_export_api",
        "optional_parameters": [
            {"default": None, "description": "Sample metadata TSV (required for ANCOM)", "name": "metadata_tsv", "type": "str"},
            {"default": ("ancom",), "description": "Differential abundance methods to attempt", "name": "methods", "type": "tuple"},
            {"default": True, "description": "Export feature table to BIOM format", "name": "export_biom", "type": "bool"},
            {"default": True, "description": "Export taxonomy to TSV format", "name": "export_tsv", "type": "bool"},
            {"default": True, "description": "Unpack .qzv visualizations", "name": "unpack_qzv", "type": "bool"},
            {"default": False, "description": "Re-run even if outputs exist", "name": "force", "type": "bool"},
        ],
        "required_parameters": [
            {"default": None, "description": "Feature table artifact (.qza)", "name": "table_qza", "type": "str"},
            {"default": None, "description": "Taxonomy artifact (.qza)", "name": "taxonomy_qza", "type": "str"},
            {"default": "./q2_stats_export_out", "description": "Directory to save outputs", "name": "output_dir", "type": "str"},
        ],
    },
    {
        "description": "Run Cactus whole-genome multiple sequence alignment, "
                       "generating hierarchical alignment (HAL) output.",
        "name": "run_cactus",
        "optional_parameters": [
            {"default": "cactus", "description": "Cactus executable name (defaults to 'cactus' in conda env)", "name": "cactus_exe", "type": "str"},
            {"default": None, "description": "Extra arguments to pass to Cactus (e.g., ['--maxCores', '8'])", "name": "extra_args", "type": "list[str]"},
        ],
        "required_parameters": [
            {"default": None, "description": "Job store directory (Cactus intermediate files)", "name": "job_store", "type": "str"},
            {"default": None, "description": "Seqfile describing genomes for alignment", "name": "seqfile", "type": "str"},
            {"default": None, "description": "Output HAL file path", "name": "out_hal", "type": "str"},
        ],
    },
    
    {
        "description": "Run Mash-based genome de-replication on a set of FASTA files. "
                    "Builds Mash sketches, computes all-vs-all distances, clusters genomes "
                    "using single-linkage under a distance cutoff (default 0.004 ≈ 99.6% identity), "
                    "and selects one representative genome per cluster (largest file by size). "
                    "Outputs include representative genome list, cluster membership, Mash distance table, "
                    "and sketch file for downstream use.",
        "name": "run_mash_deduplicate",
        "optional_parameters": [
            {"default": 21, "description": "Mash k-mer size for sketching", "name": "kmer", "type": "int"},
            {"default": 1000, "description": "Mash sketch size (-s parameter)", "name": "sketch_size", "type": "int"},
            {"default": 0.004, "description": "Maximum Mash distance for genomes to be considered redundant (≈99.6% identity)", "name": "max_distance", "type": "float"},
            {"default": 1, "description": "Number of CPU threads for Mash distance calculation", "name": "threads", "type": "int"},
            {"default": True, "description": "Keep intermediate files (temporary list, sketch, etc.)", "name": "keep_temps", "type": "bool"},
            {"default": "mash_self_dist.tsv", "description": "Filename for Mash distance table output", "name": "distance_tsv_name", "type": "str"}
        ],
        "required_parameters": [
            {"default": None, "description": "List of genome FASTA file paths to process", "name": "fasta_paths", "type": "list[str]"},
            {"default": None, "description": "Absolute path to Mash executable (do not rely on PATH)", "name": "mash_path", "type": "str"},
            {"default": "./mash_out", "description": "Output directory to save Mash results", "name": "out_dir", "type": "str"}
        ]
    },

    {
        "description": "Run ARTS (Antibiotic Resistant Target Seeker) on input FASTA to identify antibiotic resistance-related targets.",
        "name": "run_arts",
        "optional_parameters": [
            {"default": None, "description": "Extra ARTS command-line arguments", "name": "extra_args", "type": "list[str]"}
        ],
        "required_parameters": [
            {"default": None, "description": "Absolute path to ARTS executable", "name": "arts_exec", "type": "str"},
            {"default": None, "description": "Input FASTA file path", "name": "input_fasta", "type": "str"},
            {"default": "./arts_out", "description": "Output directory", "name": "output_dir", "type": "str"},
        ],
    },
    {
        "description": "Run Exonerate alignment between query and target FASTA sequences using a specified model, "
                    "saving results to an output file.",
        "name": "run_exonerate",
        "optional_parameters": [
            {"default": "exonerate", "description": "Exonerate executable name (defaults to 'exonerate' in conda env)", "name": "exonerate_exe", "type": "str"},
            {"default": None, "description": "Extra arguments to pass to Exonerate (e.g., ['--showalignment', 'yes'])", "name": "extra_args", "type": "list[str]"},
        ],
        "required_parameters": [
            {"default": None, "description": "Path to query FASTA file", "name": "query_fasta", "type": "str"},
            {"default": None, "description": "Path to target FASTA file", "name": "target_fasta", "type": "str"},
            {"default": None, "description": "Alignment model (e.g., 'protein2genome', 'dna2dna')", "name": "model", "type": "str"},
            {"default": None, "description": "Path to output result file", "name": "out_file", "type": "str"},
        ],
    },
    {
        "description": "Run FastQC on sequencing reads to assess quality metrics.",
        "name": "run_fastqc",
        "optional_parameters": [
            {"default": 1, "description": "Number of threads", "name": "threads", "type": "int"},
            {"default": None, "description": "Extra FastQC arguments", "name": "extra_args", "type": "list[str]"}
        ],
        "required_parameters": [
            {"default": None, "description": "Path to FastQC executable", "name": "fastqc_path", "type": "str"},
            {"default": None, "description": "List of input FASTQ files", "name": "inputs", "type": "list[str]"},
            {"default": "./fastqc_out", "description": "Output directory", "name": "out_dir", "type": "str"},
        ],
    },
    {
        "description": "Run Trimmomatic to trim and filter sequencing reads (single-end or paired-end).",
        "name": "run_trimmomatic",
        "optional_parameters": [
            {"default": None, "description": "Adapters FASTA file", "name": "adapters_fa", "type": "str"},
            {"default": 4, "description": "Number of threads", "name": "threads", "type": "int"},
            {"default": None, "description": "Extra Trimmomatic arguments", "name": "extra_args", "type": "list[str]"}
        ],
        "required_parameters": [
            {"default": None, "description": "Path to Java executable", "name": "java_path", "type": "str"},
            {"default": None, "description": "Path to Trimmomatic JAR file", "name": "trimmomatic_jar", "type": "str"},
            {"default": None, "description": "Mode: 'PE' or 'SE'", "name": "mode", "type": "str"},
            {"default": None, "description": "Input R1 FASTQ file", "name": "input_R1", "type": "str"},
            {"default": None, "description": "Output R1 FASTQ (or paired)", "name": "output_R1", "type": "str"},
        ],
    },
    {
        "description": "Run NCBI Datasets CLI command to download genomes, genes, proteins, or other biological data.",
        "name": "run_ncbi_datasets",
        "optional_parameters": [
            {"default": "datasets", "description": "NCBI Datasets executable (defaults to 'datasets' in conda env)", "name": "datasets_exe", "type": "str"},
        ],
        "required_parameters": [
            {"default": None, "description": "List of subcommands and arguments, e.g. ['download','genome','accession','--input-file','acc.txt','--filename','dl.zip']", "name": "subcommand", "type": "list[str]"},
            {"default": "./ncbi_datasets_out", "description": "Output directory for downloaded files", "name": "out_dir", "type": "str"},
        ],
    },
    {
        "description": "Run NCBI Dataformat CLI command to convert datasets output (JSONL) into TSV or other formats.",
        "name": "run_ncbi_dataformat",
        "optional_parameters": [
            {"default": "dataformat", "description": "NCBI Dataformat executable (defaults to 'dataformat' in conda env)", "name": "dataformat_exe", "type": "str"},
            {"default": None, "description": "Working directory to run the command in", "name": "workdir", "type": "str"},
        ],
        "required_parameters": [
            {"default": None, "description": "List of subcommands and arguments, e.g. ['tsv','genome','--fields','accession,organism-name','--inputfile','data.jsonl']", "name": "subcommand", "type": "list[str]"},
        ],
    },
    {
        "description": "Run InterProScan to annotate protein domains and functional signatures.",
        "name": "run_interproscan",
        "optional_parameters": [
            {"default": "TSV", "description": "Output formats (comma-separated)", "name": "formats", "type": "str"},
            {"default": None, "description": "Applications to run (comma-separated)", "name": "applications", "type": "list[str]"},
            {"default": None, "description": "Extra InterProScan arguments", "name": "extra_args", "type": "list[str]"}
        ],
        "required_parameters": [
            {"default": None, "description": "Path to InterProScan executable", "name": "interproscan_exec", "type": "str"},
            {"default": None, "description": "Input FASTA file", "name": "input_fasta", "type": "str"},
            {"default": "./interproscan_out", "description": "Output directory", "name": "out_dir", "type": "str"},
        ],
    },
    {
        "description": "Run Bowtie2 to align sequencing reads against a reference index, producing SAM output.",
        "name": "run_bowtie2",
        "optional_parameters": [
            {"default": 4, "description": "Number of threads", "name": "threads", "type": "int"},
            {"default": None, "description": "Extra Bowtie2 arguments", "name": "extra_args", "type": "list[str]"}
        ],
        "required_parameters": [
            {"default": None, "description": "Path to Bowtie2 executable", "name": "bowtie2_exec", "type": "str"},
            {"default": None, "description": "Index basename (prefix of .bt2 files)", "name": "index_base", "type": "str"},
            {"default": None, "description": "Read1 FASTQ file", "name": "read1", "type": "str"},
            {"default": "alignment.sam", "description": "Output SAM file", "name": "out_sam", "type": "str"},
        ],
    },
    {
        "description": "Run ANNOVAR (table_annovar.pl) to annotate variants in a VCF file against human or custom databases.",
        "name": "run_annovar",
        "optional_parameters": [
            {"default": None, "description": "Extra ANNOVAR arguments", "name": "other_args", "type": "list[str]"}
        ],
        "required_parameters": [
            {"default": None, "description": "Path to table_annovar.pl script", "name": "table_annovar_pl", "type": "str"},
            {"default": None, "description": "Path to annovar_dir (with scripts)", "name": "annovar_dir", "type": "str"},
            {"default": None, "description": "Path to humandb_dir (with annotation databases)", "name": "humandb_dir", "type": "str"},
            {"default": None, "description": "Input VCF file", "name": "input_vcf", "type": "str"},
            {"default": "hg19", "description": "Genome build version", "name": "buildver", "type": "str"},
            {"default": "anno_out", "description": "Output prefix", "name": "out_prefix", "type": "str"},
            {"default": None, "description": "Protocols list (comma-joined)", "name": "protocols", "type": "list[str]"},
            {"default": None, "description": "Operations list (comma-joined)", "name": "operations", "type": "list[str]"},
        ],
    },
]