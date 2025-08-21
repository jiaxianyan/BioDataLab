# collect_metadata

**Python Library Dependencies:**

*   **biopython** (pip/conda: `biopython`) - Used for interacting with NCBI Entrez.  Minimal version not strictly specified, but a recent version is recommended for best compatibility with NCBI's API.
*   **csv** (built-in) - No installation needed.
*   **typing** (built-in) - No installation needed.
*   **logging** (built-in) - No installation needed.

**External Programs/CLIs:**

*   None

**OS-Level Prerequisites:**

*   None

**Suggested Environment Setup:**

It's highly recommended to use a virtual environment (venv) or conda environment to manage dependencies.

*   **venv:**

    ```bash
    python3 -m venv .venv
    source .venv/bin/activate  # On Linux/macOS
    .venv\Scripts\activate  # On Windows
    pip install biopython
    ```

*   **conda:**

    ```bash
    conda create -n myenv python=3.x  # Replace 3.x with your desired Python version
    conda activate myenv
    conda install -c conda-forge biopython
    ```

**Install Commands:**

```bash
pip install biopython
```

or

```bash
conda install -c conda-forge biopython
```

**Notes:**

*   **NCBI API Key:**  Obtain an NCBI API key from the NCBI website ([https://ncbiinsights.ncbi.nlm.nih.gov/news/news-ncbi-api-keys/](https://ncbiinsights.ncbi.nlm.nih.gov/news/news-ncbi-api-keys/)) to avoid rate limiting when using the Entrez API.  Pass it to the `collect_metadata` function using the `ncbi_api_key` argument.  Without an API key, Entrez will be used without one, and you may encounter rate limiting.
*   **Email Address:**  Provide a valid email address to the `collect_metadata` function using the `email` argument. This is required by Entrez.
*   **File Paths:** Ensure that the paths to `accession_list_path`, `mammaldiet_path`, `eltontraits_path`, and `output_path` are correct and that the program has read permissions for the input files and write permissions for the output directory.
*   **Large File Handling:** The `mammaldiet.tsv` and `eltontraits.csv` files are loaded entirely into memory. If these files are very large, consider using a more memory-efficient approach, such as reading them in chunks or using a database.
*   **Error Handling:** The code includes basic error handling, but you might want to add more robust error handling, especially for network-related issues when fetching data from NCBI.  Consider adding retry logic for failed NCBI requests.
*   **Performance:** Fetching data from NCBI for each accession can be slow.  Consider implementing caching to store the results of NCBI queries to improve performance, especially if you are processing a large number of accessions.  Also, consider using asynchronous requests to NCBI to fetch data in parallel.
*   **Data Integrity:** The code assumes that the input files (`mammaldiet.tsv`, `eltontraits.csv`, `accession_list.txt`) are properly formatted. Add validation checks to ensure data integrity.
*   **Logging:** The code uses basic logging. Consider adding more detailed logging to help with debugging.
*   **Entrez Usage:** Be mindful of NCBI's Entrez API usage guidelines to avoid being blocked.


# merge_paired_fastq_with_vsearch

This function merges paired-end FASTQ files using VSEARCH. Here's a breakdown of its environment and dependencies:

**Python Library Dependencies:**

*   No explicit Python library dependencies beyond the standard library (e.g., `os`, `subprocess`, `typing`).  Therefore, no `pip install` or `conda install` commands are needed for Python libraries.

**External Programs/CLIs:**

*   **VSEARCH:**  This is the primary dependency.
    *   **Name:** `vsearch`
    *   **Minimal Version:**  While not explicitly checked, VSEARCH version 2.7.1 or later is recommended for reliable paired-end merging.  Older versions might have issues.
    *   **Install Methods:**
        *   **conda:** `conda install -c bioconda vsearch`
        *   **apt (Debian/Ubuntu):** `sudo apt-get install vsearch` (May provide an older version; check the version after installation)
        *   **brew (macOS):** `brew install vsearch`
        *   **Download and Compile:**  Download the source code from the VSEARCH website ([https://github.com/torognes/vsearch](https://github.com/torognes/vsearch)) and follow the compilation instructions. This is useful for getting the latest version or customizing the build.
    *   **Path:** The `vsearch_bin` argument *must* be the full path to the VSEARCH executable (e.g., `/usr/local/bin/vsearch`, `/opt/conda/envs/myenv/bin/vsearch`).  The function explicitly checks if the provided path exists.

**OS-Level Prerequisites:**

*   **General:**  A standard operating system (Linux, macOS, Windows with WSL) capable of running Python and VSEARCH.
*   **Linux:**  Ensure `build-essential` or equivalent packages are installed if you plan to compile VSEARCH from source.
*   **Permissions:** The user running the Python script needs execute permissions on the VSEARCH binary and read/write permissions on the input and output FASTQ files and directories.

**Suggested Environment Setup:**

It's highly recommended to use a virtual environment (venv) or Conda environment to manage dependencies.

*   **venv:**
    ```bash
    python3 -m venv .venv
    source .venv/bin/activate
    # No pip installs needed for this specific function, but you might have other dependencies in your project.
    ```

*   **conda:**
    ```bash
    conda create -n myenv python=3.9  # Or your preferred Python version
    conda activate myenv
    conda install -c bioconda vsearch
    ```

**Example Usage:**

```python
import subprocess
import os
from typing import Tuple

def merge_paired_fastq_with_vsearch(fastq1_path: str, fastq2_path: str, output_fastq_path: str, vsearch_bin: str) -> None:
    """Merges paired-end FASTQ files using VSEARCH.

    Args:
        fastq1_path: Path to the first FASTQ file (read 1).
        fastq2_path: Path to the second FASTQ file (read 2).
        output_fastq_path: Path to the output merged FASTQ file.
        vsearch_bin: Path to the VSEARCH executable.

    Returns:
        None. The merged FASTQ file is written to the specified output path.

    Raises:
        FileNotFoundError: If any of the input files or the VSEARCH executable are not found.
        subprocess.CalledProcessError: If the VSEARCH command fails.
        ValueError: If input FASTQ files are empty.

    Examples:
        merge_paired_fastq_with_vsearch(
            fastq1_path="read1.fastq",
            fastq2_path="read2.fastq",
            output_fastq_path="merged.fastq",
            vsearch_bin="/usr/local/bin/vsearch"
        )
    """

    if not os.path.isfile(fastq1_path):
        raise FileNotFoundError(f"FASTQ file not found: {fastq1_path}")
    if not os.path.isfile(fastq2_path):
        raise FileNotFoundError(f"FASTQ file not found: {fastq2_path}")
    if not os.path.isfile(vsearch_bin):
        raise FileNotFoundError(f"VSEARCH executable not found: {vsearch_bin}")

    if os.stat(fastq1_path).st_size == 0:
        raise ValueError(f"FASTQ file is empty: {fastq1_path}")
    if os.stat(fastq2_path).st_size == 0:
        raise ValueError(f"FASTQ file is empty: {fastq2_path}")

    try:
        command = [
            vsearch_bin,
            "--fastq_mergepairs", fastq1_path,
            "--reverse", fastq2_path,
            "--fastqout", output_fastq_path,
            "--fastq_minmergelen", "10" #Added to avoid vsearch error with short reads
        ]

        subprocess.run(command, check=True, capture_output=True, text=True)

    except subprocess.CalledProcessError as e:
        raise subprocess.CalledProcessError(e.returncode, e.cmd, output=e.output, stderr=e.stderr)


# Example usage (assuming VSEARCH is installed and in your PATH, or you provide the full path):
fastq1 = "read1.fastq"
fastq2 = "read2.fastq"
output = "merged.fastq"
vsearch_path = "/usr/local/bin/vsearch"  # Replace with the actual path to vsearch

# Create dummy fastq files for testing
with open(fastq1, "w") as f:
    f.write("@read1\nACGT\n+\n!!!!\n")
with open(fastq2, "w") as f:
    f.write("@read1\nTGCA\n+\n!!!!\n")

try:
    merge_paired_fastq_with_vsearch(fastq1, fastq2, output, vsearch_path)
    print(f"Successfully merged {fastq1} and {fastq2} into {output}")
except FileNotFoundError as e:
    print(f"Error: {e}")
except subprocess.CalledProcessError as e:
    print(f"VSEARCH error: {e.stderr}")
except ValueError as e:
    print(f"Value Error: {e}")
finally:
    # Clean up dummy files
    os.remove(fastq1)
    os.remove(fastq2)
    os.remove(output)
```

**Notes:**

*   **Paths:**  Be extremely careful with file paths.  Use absolute paths or paths relative to the script's execution directory to avoid confusion.
*   **Permissions:**  Ensure the script has the necessary permissions to read the input FASTQ files and write the output FASTQ file.
*   **Large Files:**  VSEARCH is generally efficient, but merging very large FASTQ files can still be memory-intensive.  Monitor memory usage if you encounter issues. Consider splitting the input files into smaller chunks if necessary (outside of this function).
*   **Error Handling:** The function includes basic error handling (file not found, VSEARCH errors, empty files).  Consider adding more robust error handling, such as logging, for production environments.
*   **Performance:** The `--fastq_minmergelen` parameter is added to the vsearch command to avoid errors with short reads. You may need to adjust this parameter based on the expected read lengths in your FASTQ files. Other VSEARCH parameters (e.g., related to quality filtering) can be added to the `command` list to further optimize the merging process.
*   **VSEARCH Output:**  VSEARCH provides detailed output to standard error.  The `subprocess.CalledProcessError` exception captures this output, which can be helpful for debugging.
*   **Alternative Merging Tools:** Other tools exist for merging paired-end reads (e.g., `flash`, `pear`). The choice of tool depends on the specific requirements of your project.


# fastq_quality_filter

This function requires the QIIME 2 environment to be set up correctly.

**Python Library Dependencies:**

*   None explicitly used in the code, but `typing` is part of the standard library.

**External Programs/CLIs:**

*   **QIIME 2:**  (Minimum version not explicitly required, but the code was written with QIIME 2 in mind, so a recent version like 2023.9 or later is recommended).
    *   Install method: QIIME 2 is typically installed using conda.  See the official QIIME 2 installation guide: [https://docs.qiime2.org/](https://docs.qiime2.org/)

**OS-Level Prerequisites:**

*   A suitable operating system for QIIME 2 (Linux or macOS are most common).
*   Conda or Miniconda installed.

**Suggested Environment Setup:**

It's highly recommended to use a dedicated conda environment for QIIME 2.

```bash
# Create a new conda environment (replace qiime2-2023.9 with your desired version)
conda create -n qiime2-2023.9 -c conda-forge -c bioconda qiime2=2023.9

# Activate the environment
conda activate qiime2-2023.9

# (Optional) Install any other Python packages you need within this environment
# pip install ...
```

**Notes:**

*   **`qiime2_bin` Path:** The `qiime2_bin` argument should point to the *directory* where the `qiime` executable is located, or directly to the `qiime` executable itself.  For example, if you installed QIIME 2 in a conda environment named `qiime2-2023.9`, and you activated that environment, the path might be something like `/home/user/miniconda3/envs/qiime2-2023.9/bin` or `/opt/qiime2-2023.9/bin`.  The code checks if the provided path is a file or a directory.
*   **Permissions:** Ensure that the user running the script has execute permissions for the QIIME 2 binary.
*   **Large File Handling:** For very large FASTQ files, consider using a more memory-efficient approach within QIIME 2 if performance is an issue.  The current implementation uses `subprocess.run` which should handle reasonably large files, but extremely large files might benefit from streaming approaches if QIIME 2 supports them.
*   **Error Handling:** The code includes basic error handling for file not found and invalid parameter values.  The `subprocess.CalledProcessError` exception provides access to the standard output and standard error from the QIIME 2 command, which can be helpful for debugging.
*   **QIIME 2 Version Compatibility:**  QIIME 2 commands and parameters can change between versions.  Ensure that the parameters used in the `fastq_quality_filter` function are compatible with the version of QIIME 2 you are using.  Consult the QIIME 2 documentation for your specific version.
*   **Alternative QIIME 2 filtering methods:** QIIME 2 offers other filtering methods. This function uses `demux filter-seqs`. Other methods might be more appropriate depending on the specific data and research question.


# denoise_with_deblur

This function denoises FASTQ files using Deblur and converts the output to a BIOM table using QIIME 2 tools. Here's a breakdown of the environment and dependencies:

**Python Library Dependencies:**

*   None explicitly required in the code, but `shutil` is used for cleanup, which is part of the Python standard library.

**External Programs/CLIs:**

*   **Deblur:** (Name: `deblur`, install method varies, often via `conda` or direct download/compilation).  The code expects the `deblur` executable to be in the `deblur_bin` directory.  You'll need to install Deblur separately.  A typical installation method is via conda: `conda install -c bioconda deblur`.  The code uses `deblur denoise-16S`.  The version of Deblur isn't explicitly checked, but the code relies on the standard command-line interface of Deblur.
*   **QIIME 2:** (Name: `qiime`, install method: `conda`). The code expects the `qiime` executable (specifically, `qiime tools export` and `qiime tools import`) to be in the `qiime2_bin` directory.  You'll need to install QIIME 2 separately. A typical installation method is via conda: `conda env create -n qiime2-2023.9 --file qiime2-2023.9-py38-linux-conda.yml` (replace `2023.9` and `py38` with your desired version and Python version).  The code uses `qiime tools export` and `qiime tools import`. The code assumes the QIIME 2 environment is activated when running the script.
*   **vsearch:** (Implicit dependency). The code creates a feature table in TSV format and imports it into QIIME 2. QIIME 2 relies on vsearch to perform sequence matching. vsearch is not explicitly called in the code, but it is a dependency of QIIME 2. Install via conda: `conda install -c bioconda vsearch`.

**OS-Level Prerequisites:**

*   Basic OS tools for file manipulation (present in most Linux/macOS environments).
*   Sufficient disk space for intermediate files (the `temp_deblur_output` directory).

**Suggested Environment Setup:**

1.  **Conda Environment:**  This is highly recommended due to the complexity of QIIME 2 and Deblur dependencies.

    ```bash
    conda create -n deblur_qiime python=3.8  # Or your preferred Python version
    conda activate deblur_qiime
    conda install -c bioconda deblur qiime2 # Installs latest versions.  Consider specifying versions.
    conda install -c bioconda vsearch
    ```

    *   Replace `python=3.8` with the Python version compatible with your QIIME 2 version.
    *   Consider specifying the exact versions of `deblur` and `qiime2` to ensure reproducibility.  For example, `conda install -c bioconda deblur=1.1.1 qiime2=2023.9`.

2.  **Alternative: venv (less recommended):**  While possible, managing the dependencies for QIIME 2 and Deblur with `venv` and `pip` is significantly more complex and error-prone.  It's strongly advised to use Conda.

**Example Install Commands (within the Conda environment):**

```bash
# After creating and activating the conda environment (see above)
# Deblur and QIIME 2 are installed as part of the environment creation.
# vsearch is installed as part of the environment creation.
```

**Notes:**

*   **Paths:** The function relies on the correct paths to the `qiime2_bin` and `deblur_bin` directories.  Ensure these are set correctly.  Absolute paths are recommended for clarity and to avoid issues with relative paths.
*   **Permissions:**  The user running the script needs read permissions for the input FASTQ file and execute permissions for the `deblur` and `qiime` executables.  Write permissions are needed for the `output_biom` file and the temporary directory.
*   **Large File Handling:**  For very large FASTQ files, consider using a more memory-efficient method for reading and processing the sequences within the temporary directory creation. The current implementation reads the entire FASTA file into memory.
*   **Performance:**
    *   The `--threads` option for Deblur is currently hardcoded to `1`.  Consider making this configurable to leverage multi-core processors for faster denoising.  However, be mindful of memory usage when increasing the number of threads.
    *   Ensure that the input FASTQ file is efficiently compressed (e.g., using `gzip`).
*   **Error Handling:** The code includes basic error handling for file not found and subprocess errors. Consider adding more specific error handling and logging for debugging purposes.
*   **Deblur Databases:** The code assumes that `core_alignment.fasta` and `negative_control_sequences.fasta` are located in the `deblur_bin` directory. Ensure these files are present and are the correct versions for your Deblur installation.
*   **QIIME 2 Version:** The code is written assuming a relatively recent version of QIIME 2. Older versions might have different command-line arguments or API behavior.
*   **Temporary Directory:** The temporary directory `temp_deblur_output` is created in the current working directory. Consider making this configurable and using a more robust method for creating temporary directories (e.g., `tempfile.mkdtemp()`).
*   **Sequence Type:** The code imports the feature table as `FeatureData[Sequence]`. If you need to import it as a different type (e.g., `FeatureData[AlignedSequence]`), you'll need to modify the `--type` argument in the `qiime tools import` command.


# filter_biom_table

**Python Library Dependencies:**

*   `biom-format`: Required for reading, filtering, and writing BIOM tables. Install with `pip install biom-format`.  A specific minimum version isn't strictly required, but it's generally good practice to use a relatively recent version (e.g., 2.1.7 or later) to benefit from bug fixes and performance improvements.

**External Programs/CLIs:**

*   None. This function relies solely on the `biom-format` Python library.

**OS-Level Prerequisites:**

*   None. This function is platform-independent as long as Python and the `biom-format` library are installed.

**Suggested Environment Setup:**

It's highly recommended to use a virtual environment (venv) or conda environment to manage dependencies.

*   **venv:**

    ```bash
    python3 -m venv .venv
    source .venv/bin/activate  # On Linux/macOS
    # .venv\Scripts\activate  # On Windows
    pip install biom-format
    ```

*   **conda:**

    ```bash
    conda create -n biom_env python=3.x  # Replace 3.x with your desired Python version
    conda activate biom_env
    conda install -c conda-forge biom-format
    ```

**Notes:**

*   **Paths:** The function expects valid file paths for both the input BIOM table and the output file. Ensure the input file exists and the output directory has write permissions.
*   **Permissions:** The user running the script needs read permissions for the input BIOM table and write permissions for the output directory.
*   **Large-File Handling:** The `biom-format` library is designed to handle large BIOM tables efficiently. However, very large tables might still require significant memory. If you encounter memory issues, consider processing the table in chunks or using a more memory-efficient approach if available within the `biom-format` library (check the library's documentation for advanced usage).
*   **Performance Tips:** For very large BIOM tables, consider using optimized versions of NumPy and SciPy, which `biom-format` relies on. These can be installed via conda or pip with appropriate flags for optimized builds.
*   **Error Handling:** The function includes basic error handling for file existence, invalid input types, and BIOM table processing errors.  More specific error handling might be needed depending on the application.
*   **Output Format:** The output BIOM table is saved in HDF5 format, which is the standard format for BIOM tables.


# perform_diversity_analysis

**Environment Needs:**

*   **Python Library Dependencies:**
    *   No explicit Python library dependencies beyond the standard library (e.g., `os`, `subprocess`, `typing`).

*   **External Programs/CLIs:**
    *   **QIIME 2:** (Minimum version not explicitly required, but the code assumes a QIIME 2 installation with the `qiime` binary accessible).
        *   Install Method: QIIME 2 is typically installed using conda.  See the official QIIME 2 installation guide: [https://docs.qiime2.org/](https://docs.qiime2.org/)

*   **OS-Level Prerequisites:**
    *   A functional operating system (Linux, macOS, or Windows via WSL) capable of running QIIME 2.
    *   Sufficient disk space for the ASV table, intermediate files, and output files.

**Suggested Environment Setup:**

It is highly recommended to use a dedicated conda environment for QIIME 2 and this script. This isolates the QIIME 2 installation and avoids conflicts with other Python packages.

1.  **Create a conda environment:**

    ```bash
    conda create -n qiime2-env -c conda-forge python=3.8  # Or a compatible Python version
    conda activate qiime2-env
    ```

2.  **Install QIIME 2:**

    Follow the official QIIME 2 installation instructions for your operating system.  A typical installation might look like this (but refer to the QIIME 2 documentation for the most up-to-date instructions):

    ```bash
    conda install -c conda-forge -c bioconda qiime2
    ```

3.  **Install any other necessary Python packages (if any):**

    If you add any other python dependencies to the script, install them using pip:

    ```bash
    pip install <package_name>
    ```

**Notes:**

*   **`qiime2_bin` Path:** The `qiime2_bin` argument is crucial.  It *must* point to the correct path of the `qiime` executable within your QIIME 2 environment.  Verify this path after installing QIIME 2.  A common location is `/opt/qiime2-<version>/bin/qiime` or within your conda environment's `bin` directory (e.g., `~/miniconda3/envs/qiime2-env/bin/qiime`).

*   **Permissions:** Ensure that the user running the script has execute permissions for the `qiime2_bin` and write permissions for the `output_dir`.

*   **Large File Handling:**  ASV tables can be very large.  Ensure you have enough RAM to handle the data, especially during rarefaction.  Consider using a server or a machine with substantial memory if you encounter memory errors.

*   **Metadata File:** The code uses the ASV table as the metadata file for alpha and beta group significance. This is likely incorrect. A separate metadata file is usually required, containing sample IDs and associated metadata (e.g., treatment groups, environmental variables).  You'll need to create a proper metadata file and update the `perform_diversity_analysis` function to use it. The metadata file should be a tab-separated file.

*   **Error Handling:** The code includes basic error handling for file not found and invalid rarefaction depth.  It also captures the output and stderr from QIIME 2 commands, which is helpful for debugging.

*   **QIIME 2 Version Compatibility:** QIIME 2 versions can introduce changes that might affect command syntax or output formats.  Test the script thoroughly with your specific QIIME 2 version.

*   **Performance:** Rarefaction can be computationally intensive.  Consider using a smaller rarefaction depth for initial testing or if performance is a concern.

*   **Output Files:** The function creates several output files in the `output_dir`:
    *   `rarefied_table.biom`: The rarefied ASV table.
    *   `alpha_vector.qza`: Alpha diversity results.
    *   `beta_diversity.qza`: Beta diversity distance matrix.
    *   `alpha_visualization.qzv`: Visualization of alpha group significance.
    *   `beta_visualization.qzv`: Visualization of beta group significance.
    The `.qza` files are QIIME 2 artifacts, and the `.qzv` files are QIIME 2 visualizations that can be viewed using `qiime tools view`.


# construct_phylogenetic_tree

**Environment Needs:**

*   **Python Library Dependencies:**
    *   None explicitly required beyond the Python standard library.

*   **External Programs/CLIs:**
    *   **MAFFT:** (Minimal version not strictly defined, but recent versions are recommended for accuracy and speed. v7 is widely used).
        *   Install methods:
            *   **conda:** `conda install -c bioconda mafft`
            *   **apt (Debian/Ubuntu):** `sudo apt-get install mafft`
            *   **brew (macOS):** `brew install mafft`
            *   From source: Download from MAFFT website and compile.
    *   **FastTree:** (Minimal version not strictly defined, but recent versions are recommended for accuracy and speed. v2 is common).
        *   Install methods:
            *   **conda:** `conda install -c bioconda fasttree`
            *   **apt (Debian/Ubuntu):** `sudo apt-get install fasttree`
            *   **brew (macOS):** `brew install fasttree`
            *   From source: Download from FastTree website and compile.

*   **OS-level prerequisites:**
    *   A functional operating system (Linux, macOS, Windows via WSL).
    *   A C++ compiler might be needed if installing MAFFT or FastTree from source.

*   **Suggested environment setup:**

    It's highly recommended to use a virtual environment (venv) or Conda environment to manage dependencies.

    *   **venv:**
        ```bash
        python3 -m venv .venv
        source .venv/bin/activate  # On Linux/macOS
        # .venv\Scripts\activate  # On Windows
        # No pip installs needed, as there are no explicit python dependencies.
        ```

    *   **conda:**
        ```bash
        conda create -n phylo_env python=3.8  # Or any desired Python version
        conda activate phylo_env
        conda install -c bioconda mafft fasttree
        ```

*   **Notes:**

    *   **Paths:** The `mafft_bin` and `fasttree_bin` arguments are crucial.  Ensure these paths are correct and point to the actual executables.  Use `which mafft` and `which FastTree` (or `where mafft` and `where FastTree` on Windows) to find the full paths after installation.
    *   **Permissions:** The user running the script needs execute permissions for the MAFFT and FastTree binaries.  They also need write permissions to the `output_dir`.
    *   **Large-file handling:** For very large datasets, MAFFT can be memory-intensive. Consider using the `--thread` option with MAFFT to utilize multiple cores, if available.  Also, monitor memory usage.  FastTree is generally less memory-intensive.
    *   **Performance tips:**
        *   Use the `--thread` option in MAFFT if you have multiple CPU cores.  The `--auto` option in MAFFT usually selects a reasonable alignment strategy, but for very large datasets, you might need to experiment with different alignment algorithms (e.g., `--localpair`).
        *   For extremely large datasets, consider using a more scalable phylogenetic tree building method than FastTree, such as RAxML-NG or IQ-TREE, but these are more complex to set up and use.
    *   **Error Handling:** The code includes basic error handling for missing executables, empty sequence lists, and failed subprocess calls.  Consider adding more robust error handling, such as checking the validity of the ASV sequences (e.g., ensuring they only contain valid nucleotide characters).
    *   **Input Validation:** The code checks for empty ASV sequences, but it doesn't validate the format of the sequences themselves. Consider adding validation to ensure that the sequences are valid FASTA format and contain only valid nucleotide characters (A, T, G, C, N).
    *   **Temporary Files:** The code creates a temporary FASTA file (`asv_sequences.fasta`). Consider using the `tempfile` module to create a truly temporary file that is automatically deleted when the script finishes. This can help prevent clutter in the output directory.


# annotate_asv_with_qiime2

This function annotates ASV sequences using QIIME2's `classify-consensus-vsearch` method against a reference database (typically EzBioCloud). Here's a breakdown of the environment and dependencies:

**Python Library Dependencies:**

*   None explicitly required beyond the Python standard library (e.g., `os`, `subprocess`).  No `pip install` or `conda install` commands are needed for *this specific function* beyond the QIIME2 installation itself.

**External Programs/CLIs:**

*   **QIIME2:**  This is the primary dependency.
    *   Name: `qiime`
    *   Minimal Version:  While the code doesn't enforce a specific version, it's designed with QIIME2 versions 2023.9 or later in mind.  Older versions might have different command structures or parameter names.
    *   Install Method: QIIME2 is typically installed using conda.  See the official QIIME2 installation guide: [https://docs.qiime2.org/](https://docs.qiime2.org/)
    *   Crucially, the `qiime2_bin` argument *must* point to the full path of the `qiime` executable within your QIIME2 environment (e.g., `/opt/qiime2-2023.9/bin/qiime` or `/home/user/miniconda3/envs/qiime2-2023.9/bin/qiime`).

*   **VSEARCH:** QIIME2's `classify-consensus-vsearch` method relies on VSEARCH for sequence alignment. VSEARCH is installed as part of the QIIME2 installation process, so you don't need to install it separately *if* you followed the QIIME2 installation instructions.

**OS-Level Prerequisites:**

*   A working operating system (Linux, macOS are most common for QIIME2).
*   Sufficient disk space for the QIIME2 installation, the reference database, and intermediate files.  The EzBioCloud database can be quite large.
*   Sufficient RAM.  Sequence alignment can be memory-intensive, especially with large datasets.

**Suggested Environment Setup:**

1.  **Conda Environment:**  QIIME2 *strongly* recommends using a dedicated conda environment.

2.  **Create and Activate Environment:**

    ```bash
    conda create -n qiime2-2023.9 --yes python=3.8  # Or your preferred Python version
    conda activate qiime2-2023.9
    ```

3.  **Install QIIME2:**  Follow the official QIIME2 installation instructions for your operating system.  A typical installation command looks like this:

    ```bash
    conda install -c conda-forge qiime2-2023.9
    ```

    (Replace `2023.9` with the desired QIIME2 version).

4.  **Install q2cli (if not already installed):**
    ```bash
    conda install -c conda-forge q2cli
    ```

**Notes:**

*   **Paths:**  Double-check that the paths to `asv_fasta`, `qiime2_bin`, and `ezbio_db` are correct and accessible.  Use absolute paths to avoid confusion.  The function includes checks for file existence, but it's still good practice to verify paths.
*   **Permissions:**  Ensure that the user running the script has read permissions for the input files (`asv_fasta`, `ezbio_db`) and write permissions for the output directory where `output_tsv` will be created.
*   **EzBioCloud Database:** The `ezbio_db` argument can be either a directory containing `taxonomy.qza` and `sequences.qza` *or* the direct path to the `taxonomy.qza` file. The code attempts to handle both cases.  However, it's best practice to ensure that both `taxonomy.qza` and `sequences.qza` are present and that the path to `ezbio_db` is correct.  Download the EzBioCloud database from the official source and extract it appropriately.
*   **Large File Handling:**  For very large ASV FASTA files, consider using a more memory-efficient FASTA parser if the script becomes slow.  However, the QIIME2 command itself is likely to be the bottleneck.
*   **Performance Tips:**
    *   The `--threads 1` argument is included to ensure deterministic behavior.  Removing this argument will allow QIIME2 to use multiple threads, potentially speeding up the process, but the results might vary slightly between runs.  If reproducibility is critical, leave `--threads 1` in place.
    *   Ensure that VSEARCH is properly installed and configured within the QIIME2 environment.
    *   Consider using a smaller subset of the ASV sequences for testing before running the full dataset.
*   **Error Handling:** The function includes basic error handling (file existence checks, identity threshold validation, and `subprocess.CalledProcessError` handling).  However, you might want to add more robust error handling, such as logging, to capture more detailed information about any failures.
*   **Database Format:** The code assumes the EzBioCloud database is in QIIME2 archive (`.qza`) format. If it's in a different format, you'll need to convert it to `.qza` using QIIME2's import functionality *before* using this function.
