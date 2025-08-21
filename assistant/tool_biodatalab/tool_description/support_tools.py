description = [

    {
        "description": "Executes the provided Python command in the virtual and isolated notebook environment and returns the output.",
        "name": "run_python_repl_virtual",
        "optional_parameters": [],
        "required_parameters": [
            {
                "default": None,
                "description": "Python command to execute in the notebook environment",
                "name": "command",
                "type": "str",
            }
        ],
    },
    {
        "description": "Executes the provided Python command in a real environment and returns the output",
        "name": "run_python_repl_real",
        "optional_parameters": [],
        "required_parameters": [
            {
                "default": None,
                "description": "Python command to execute in the notebook environment",
                "name": "command",
                "type": "str",
            }
        ],
    },
    {
        "description": "Read the source code of a function from any module path.",
        "name": "read_function_source_code",
        "optional_parameters": [],
        "required_parameters": [
            {
                "default": None,
                "description": "Fully qualified function name "
                "(e.g., "
                "'bioagentos.tool.support_tools.write_python_code')",
                "name": "function_name",
                "type": "str",
            }
        ],
    },
    
    {
        "description": "Downloads genome FASTA sequences from the NCBI Assembly database for a given taxonomic name, extracts the files, and generates an accession list.",
        "name": "download_ncbi_genomes_by_taxon",
        "optional_parameters": [
            {
                "default": "work/ncbi_download",
                "description": "Output directory to store the downloaded data and metadata",
                "name": "out_dir",
                "type": "str"
            },
            {
                "default": "datasets",
                "description": "Path to the NCBI 'datasets' CLI executable (default assumes it is available in the system PATH)",
                "name": "datasets_path",
                "type": "str"
            }
        ],
        "required_parameters": [
            {
                "default": None,
                "description": "Taxonomic name to query, e.g. 'cyanobacteria'",
                "name": "taxon",
                "type": "str"
            }
        ]
    },

]
