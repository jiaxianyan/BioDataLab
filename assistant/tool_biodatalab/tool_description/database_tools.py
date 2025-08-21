desciption = [    
    {
        "description": "Searches PubMed using a query string and date range, retrieves abstracts for matching PMIDs, and saves them to a TSV file.",
        "name": "pubmed_fetch_abstracts",
        "optional_parameters": [
            {
                "name": "out_tsv",
                "type": "str",
                "default": None,
                "description": "Optional path to save the TSV file containing PMIDs, titles, and abstracts. If None, results will be saved under benchmark/results/3 by default."
            }
        ],
        "required_parameters": [
            {
                "name": "query",
                "type": "str",
                "description": "PubMed search query string, e.g., 'SARS-CoV-2 AND immune escape'."
            },
            {
                "name": "date_from",
                "type": "str",
                "description": "Start date for the search, format 'YYYY/MM/DD'."
            },
            {
                "name": "date_to",
                "type": "str",
                "description": "End date for the search, format 'YYYY/MM/DD'."
            }
        ]
    },

    {
        "description": "Download SRA runs as FASTQ files using SRA Toolkit, based on query_sra results. "
                       "This will run 'prefetch' to retrieve the SRA run and 'fasterq-dump' to convert to FASTQ files.",
        "name": "download_sra_runs",
        "optional_parameters": [
            {"default": "./operation_env/tool_lake/sra_toolkit", "description": "Path to SRA Toolkit binaries (must contain prefetch and fasterq-dump)", "name": "sra_toolkit_dir", "type": "str"},
            {"default": "./sra_fastq", "description": "Directory to save downloaded FASTQ files", "name": "outdir", "type": "str"}
        ],
        "required_parameters": [
            {"default": None, "description": "Result dictionary from query_sra() containing 'formatted_results' with SRA Run IDs", "name": "query_result", "type": "dict"}
        ],
    },
    {
        "description": "Download GEO data (GSE/GSM) from NCBI FTP based on query_geo results. "
                       "Supports downloading supplementary files, SOFT, or MINiML formats from the GEO FTP hierarchy.",
        "name": "download_geo",
        "optional_parameters": [
            {"default": "./geo_data", "description": "Directory to save downloaded files", "name": "outdir", "type": "str"},
            {"default": "suppl", "description": "Which GEO FTP subdirectory to download from: 'suppl', 'soft', or 'miniml'", "name": "file_type", "type": "str"}
        ],
        "required_parameters": [
            {"default": None, "description": "Result dictionary from query_geo() containing 'formatted_results' with GEO IDs (GSE, GSM)", "name": "query_result", "type": "dict"}
        ],
    },
]