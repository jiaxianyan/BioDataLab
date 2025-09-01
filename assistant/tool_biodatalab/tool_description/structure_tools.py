description = [
    {
        "description": "Import FASTQ data via a QIIME 2 manifest, optionally join paired-end reads, "
                       "trim adapters, perform quality-score filtering, and generate demux summary.",
        "name": "q2_ingest_and_qc",
        "optional_parameters": [
            {"default": None, "description": "Sample metadata TSV for later visualizations/analyses", "name": "sample_metadata_tsv", "type": "str"},
            {"default": True, "description": "Whether reads are paired-end", "name": "paired_end", "type": "bool"},
            {"default": False, "description": "If paired_end, join pairs using VSEARCH", "name": "join_paired", "type": "bool"},
            {"default": None, "description": "Extra cutadapt parameters string (e.g. '--p-front-f ACTG --p-front-r ACTG')", "name": "cutadapt_params", "type": "str"},
            {"default": False, "description": "Enable quality-filter q-score", "name": "qscore_filter", "type": "bool"},
            {"default": None, "description": "Forward truncation length hint (recorded in provenance)", "name": "trunc_len_f", "type": "int"},
            {"default": None, "description": "Reverse truncation length hint (recorded in provenance)", "name": "trunc_len_r", "type": "int"},
            {"default": 33, "description": "Phred offset used in manifest format name (33 or 64)", "name": "phred_offset", "type": "int"},
            {"default": False, "description": "Re-run even if outputs exist", "name": "force", "type": "bool"},
        ],
        "required_parameters": [
            {"default": None, "description": "Path to QIIME 2 manifest CSV (single or paired)", "name": "manifest_csv", "type": "str"},
            {"default": "./q2_ingest_and_qc_out", "description": "Directory to save outputs", "name": "output_dir", "type": "str"},
            {"default": None, "description": "Absolute path to the 'qiime' executable", "name": "qiime_bin", "type": "str"},
        ],
    },
    {
        "description": "Denoise reads (Deblur or DADA2), generate feature table and representative sequences, "
                       "and optionally filter samples by minimal reads.",
        "name": "q2_denoise_and_feature_table",
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
            {"default": None, "description": "Absolute path to the 'qiime' executable", "name": "qiime_bin", "type": "str"},
        ],
    },
    {
        "description": "Assign taxonomy using a pretrained classifier or reference reads+taxonomy (sklearn or vsearch), "
                       "then produce taxa barplot.",
        "name": "q2_taxonomy",
        "optional_parameters": [
            {"default": None, "description": "Sample metadata TSV for taxa barplot grouping", "name": "metadata_tsv", "type": "str"},
            {"default": None, "description": "Pretrained classifier artifact (.qza)", "name": "pretrained_classifier_qza", "type": "str"},
            {"default": None, "description": "Reference reads artifact (.qza) for training or vsearch", "name": "ref_seqs_qza", "type": "str"},
            {"default": None, "description": "Reference taxonomy artifact (.qza) for training or vsearch", "name": "ref_taxonomy_qza", "type": "str"},
            {"default": "sklearn", "description": "Classification method: 'sklearn' or 'vsearch'", "name": "method", "type": "str"},
            {"default": 0.97, "description": "Identity threshold for vsearch consensus", "name": "perc_identity", "type": "float"},
            {"default": False, "description": "Re-run even if outputs exist", "name": "force", "type": "bool"},
        ],
        "required_parameters": [
            {"default": None, "description": "Feature table artifact (.qza)", "name": "table_qza", "type": "str"},
            {"default": None, "description": "Representative sequences artifact (.qza)", "name": "rep_seqs_qza", "type": "str"},
            {"default": "./q2_taxonomy_out", "description": "Directory to save outputs", "name": "output_dir", "type": "str"},
            {"default": None, "description": "Absolute path to the 'qiime' executable", "name": "qiime_bin", "type": "str"},
        ],
    },
    {
        "description": "Build multiple sequence alignment and phylogenetic tree, then run core-metrics-phylogenetic "
                       "to compute alpha/beta diversity, PCoA, and Emperor visualizations.",
        "name": "q2_phylogeny_and_core_diversity",
        "optional_parameters": [
            {"default": 1000, "description": "Rarefaction depth for core metrics", "name": "sampling_depth", "type": "int"},
            {"default": True, "description": "Use midpoint root instead of explicit rooting", "name": "use_midpoint_root", "type": "bool"},
            {"default": False, "description": "Re-run even if outputs exist", "name": "force", "type": "bool"},
        ],
        "required_parameters": [
            {"default": None, "description": "Feature table artifact (.qza)", "name": "table_qza", "type": "str"},
            {"default": None, "description": "Representative sequences artifact (.qza)", "name": "rep_seqs_qza", "type": "str"},
            {"default": "./q2_diversity_out", "description": "Directory to save outputs", "name": "output_dir", "type": "str"},
            {"default": None, "description": "Absolute path to the 'qiime' executable", "name": "qiime_bin", "type": "str"},
        ],
    },
    {
        "description": "Run common differential abundance placeholders (ANCOM/ALDEx2/Songbird notes) and export tables. "
                       "Exports BIOM and taxonomy TSV; records guidance for dataset-specific steps.",
        "name": "q2_stats_and_export",
        "optional_parameters": [
            {"default": None, "description": "Sample metadata TSV (required for most differential tests)", "name": "metadata_tsv", "type": "str"},
            {"default": ["ancom"], "description": "Requested methods: any of ['ancom','aldex2','songbird']", "name": "methods", "type": "list[str]"},
            {"default": True, "description": "Export feature table as BIOM", "name": "export_biom", "type": "bool"},
            {"default": True, "description": "Export taxonomy as TSV", "name": "export_tsv", "type": "bool"},
            {"default": True, "description": "Record instructions to unpack .qzv visualizations", "name": "unpack_qzv", "type": "bool"},
            {"default": False, "description": "Re-run even if outputs exist", "name": "force", "type": "bool"},
        ],
        "required_parameters": [
            {"default": None, "description": "Feature table artifact (.qza)", "name": "table_qza", "type": "str"},
            {"default": None, "description": "Taxonomy artifact (.qza)", "name": "taxonomy_qza", "type": "str"},
            {"default": "./q2_export_out", "description": "Directory to save outputs", "name": "output_dir", "type": "str"},
            {"default": None, "description": "Absolute path to the 'qiime' executable", "name": "qiime_bin", "type": "str"},
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
        "description": "Run Exonerate to align query sequences to a target database with a specified model.",
        "name": "run_exonerate",
        "optional_parameters": [
            {"default": None, "description": "Extra Exonerate arguments", "name": "extra_args", "type": "list[str]"}
        ],
        "required_parameters": [
            {"default": None, "description": "Path to Exonerate executable", "name": "exonerate_path", "type": "str"},
            {"default": None, "description": "Query FASTA file", "name": "query_fasta", "type": "str"},
            {"default": None, "description": "Target FASTA file", "name": "target_fasta", "type": "str"},
            {"default": "protein2genome", "description": "Exonerate model (e.g. protein2genome)", "name": "model", "type": "str"},
            {"default": "exonerate_out.txt", "description": "Output alignment file", "name": "out_file", "type": "str"},
        ],
    },
    {
        "description": "Run Cactus for multiple genome alignment, producing a HAL alignment file.",
        "name": "run_cactus",
        "optional_parameters": [
            {"default": None, "description": "Extra Cactus arguments", "name": "extra_args", "type": "list[str]"}
        ],
        "required_parameters": [
            {"default": None, "description": "Path to Cactus executable", "name": "cactus_path", "type": "str"},
            {"default": None, "description": "Job store identifier (e.g. ./jobStore)", "name": "job_store", "type": "str"},
            {"default": None, "description": "Seqfile describing genomes and relationships", "name": "seqfile", "type": "str"},
            {"default": "out.hal", "description": "Output HAL alignment file", "name": "out_hal", "type": "str"},
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
        "description": "Use NCBI 'datasets' CLI to download genomic data or metadata.",
        "name": "run_ncbi_datasets",
        "optional_parameters": [],
        "required_parameters": [
            {"default": None, "description": "Path to 'datasets' executable", "name": "datasets_exec", "type": "str"},
            {"default": None, "description": "Subcommand list (e.g. ['download','genome','accession','--input-file','acc.txt'])", "name": "subcommand", "type": "list[str]"},
            {"default": "./ncbi_datasets_out", "description": "Output directory", "name": "out_dir", "type": "str"},
        ],
    },
    {
        "description": "Use NCBI 'dataformat' CLI to convert datasets JSONL to TSV or other formats.",
        "name": "run_ncbi_dataformat",
        "optional_parameters": [],
        "required_parameters": [
            {"default": None, "description": "Path to 'dataformat' executable", "name": "dataformat_exec", "type": "str"},
            {"default": None, "description": "Subcommand list (e.g. ['tsv','genome','--fields','accession'])", "name": "subcommand", "type": "list[str]"},
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