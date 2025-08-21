

import csv
import logging
from typing import List, Dict, Optional, Tuple
from Bio import Entrez


def collect_metadata(
    accession_list_path: str,
    mammaldiet_path: str,
    eltontraits_path: str,
    output_path: str,
    ncbi_api_key: Optional[str] = None,
    email: str = "your_email@example.com",
) -> str:
    """Collects metadata for accessions from various sources and outputs a TSV file.

    Args:
        accession_list_path: Path to a text file containing a list of accessions (one per line).
        mammaldiet_path: Path to the MammalDIET TSV file.
        eltontraits_path: Path to the EltonTraits CSV file.
        output_path: Path to the output metadata TSV file.
        ncbi_api_key: Optional NCBI API key. If None, Entrez will be used without a key.
        email: Email address to use with Entrez.

    Returns:
        Path to the output metadata TSV file.

    Raises:
        FileNotFoundError: If any of the input files are not found.
        ValueError: If there are issues with the input data or API key.
        RuntimeError: If there are issues fetching data from NCBI.

    Examples:
        collect_metadata(
            "accessions.txt",
            "mammaldiet.tsv",
            "eltontraits.csv",
            "metadata.tsv",
            ncbi_api_key="YOUR_API_KEY",
            email="your_email@example.com",
        )
    """

    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    # Preflight checks
    for path in [accession_list_path, mammaldiet_path, eltontraits_path]:
        if not os.path.exists(path):
            raise FileNotFoundError(f"Input file not found: {path}")

    if not email:
        raise ValueError("Email address is required for Entrez.")

    Entrez.email = email
    if ncbi_api_key:
        Entrez.api_key = ncbi_api_key

    # Read accessions from file
    try:
        with open(accession_list_path, "r") as f:
            accessions = [line.strip() for line in f]
    except Exception as e:
        raise RuntimeError(f"Error reading accession list: {e}")

    # Load diet data
    mammaldiet_data = _load_mammaldiet(mammaldiet_path)
    eltontraits_data = _load_eltontraits(eltontraits_path)

    # Prepare output data
    metadata = []
    for accession in accessions:
        try:
            accession_data = _fetch_accession_data(accession)
            host_scientific_name = accession_data.get("HostScientificName")
            project = accession_data.get("Project")
            country = accession_data.get("Country")
            pmid = accession_data.get("PMID")

            diet_info = _merge_diet_info(host_scientific_name, mammaldiet_data, eltontraits_data)
            diet_type = diet_info.get("DietType")
            trophic_level = diet_info.get("TrophicLevel")

            metadata.append({
                "Accession": accession,
                "HostScientificName": host_scientific_name,
                "Project": project,
                "Country": country,
                "PMID": pmid,
                "DietType": diet_type,
                "TrophicLevel": trophic_level,
            })

        except Exception as e:
            logging.error(f"Error processing accession {accession}: {e}")

    # Write output to TSV
    try:
        with open(output_path, "w", newline="") as tsvfile:
            fieldnames = ["Accession", "HostScientificName", "Project", "Country", "PMID", "DietType", "TrophicLevel"]
            writer = csv.DictWriter(tsvfile, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(metadata)
    except Exception as e:
        raise RuntimeError(f"Error writing output file: {e}")

    return output_path


import os

def _load_mammaldiet(mammaldiet_path: str) -> Dict[str, Dict[str, str]]:
    """Loads data from the MammalDIET TSV file.

    Args:
        mammaldiet_path: Path to the MammalDIET TSV file.

    Returns:
        A dictionary mapping scientific names to diet information.
    """
    data = {}
    try:
        with open(mammaldiet_path, "r", newline="") as tsvfile:
            reader = csv.DictReader(tsvfile, delimiter="\t")
            for row in reader:
                scientific_name = row["ScientificName"]
                data[scientific_name] = {"DietType": row["DietType"]}
    except Exception as e:
        logging.warning(f"Error loading MammalDIET data: {e}")
    return data


def _load_eltontraits(eltontraits_path: str) -> Dict[str, Dict[str, str]]:
    """Loads data from the EltonTraits CSV file.

    Args:
        eltontraits_path: Path to the EltonTraits CSV file.

    Returns:
        A dictionary mapping scientific names to trophic level information.
    """
    data = {}
    try:
        with open(eltontraits_path, "r", newline="") as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                scientific_name = row["ScientificName"]
                data[scientific_name] = {"TrophicLevel": row["TrophicLevel"]}
    except Exception as e:
        logging.warning(f"Error loading EltonTraits data: {e}")
    return data


def _fetch_accession_data(accession: str) -> Dict[str, str]:
    """Fetches accession data from NCBI using the Entrez API.

    Args:
        accession: The accession ID.

    Returns:
        A dictionary containing accession metadata.

    Raises:
        RuntimeError: If there are issues fetching data from NCBI.
    """
    try:
        handle = Entrez.esearch(db="nucleotide", term=accession, retmax="1")
        record = Entrez.read(handle)
        handle.close()
        if not record["IdList"]:
            logging.warning(f"No record found for accession {accession}")
            return {}
        gi = record["IdList"][0]

        handle = Entrez.efetch(db="nucleotide", id=gi, rettype="gb", retmode="xml")
        records = Entrez.read(handle)
        handle.close()

        if not records:
            logging.warning(f"No GenBank record found for accession {accession}")
            return {}

        record = records[0]
        features = record.get("GBSeq_feature-table", [])
        host = None
        country = None
        project = None
        pmid = None

        for feature in features:
            if feature.get("GBFeature_key") == "source":
                for qualifier in feature.get("GBFeature_quals", []):
                    if qualifier.get("GBQualifier_name") == "host":
                        host = qualifier.get("GBQualifier_value")
                    elif qualifier.get("GBQualifier_name") == "country":
                        country = qualifier.get("GBQualifier_value")

            elif feature.get("GBFeature_key") == "Pub":
                for qualifier in feature.get("GBFeature_quals", []):
                    if qualifier.get("GBQualifier_name") == "pubmed":
                        pmid = qualifier.get("GBQualifier_value")

        project = record.get("GBSeq_project", [None])[0]

        return {
            "HostScientificName": host,
            "Project": project,
            "Country": country,
            "PMID": pmid,
        }

    except Exception as e:
        raise RuntimeError(f"Error fetching data from NCBI for accession {accession}: {e}")


def _merge_diet_info(
    host_scientific_name: Optional[str],
    mammaldiet_data: Dict[str, Dict[str, str]],
    eltontraits_data: Dict[str, Dict[str, str]],
) -> Dict[str, str]:
    """Merges diet information from MammalDIET and EltonTraits data.

    Args:
        host_scientific_name: The scientific name of the host.
        mammaldiet_data: Data loaded from MammalDIET.
        eltontraits_data: Data loaded from EltonTraits.

    Returns:
        A dictionary containing merged diet information.
    """
    if not host_scientific_name:
        return {}

    diet_info = {}
    if host_scientific_name in mammaldiet_data:
        diet_info["DietType"] = mammaldiet_data[host_scientific_name].get("DietType")
    if host_scientific_name in eltontraits_data:
        diet_info["TrophicLevel"] = eltontraits_data[host_scientific_name].get("TrophicLevel")

    return diet_info

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

import subprocess
import os
from typing import Optional


def fastq_quality_filter(
    merged_fastq: str,
    quality_threshold: int,
    min_percent_length: float,
    qiime2_bin: str,
    filtered_fastq: str,
) -> None:
    """Filters a merged FASTQ file based on quality scores using QIIME 2.

    This function uses the `qiime demux filter-seqs` command to remove low-quality reads
    from a FASTQ file. It filters reads based on a quality score threshold, minimum
    percentage of high-quality bases, and removes reads containing Ns.

    Args:
        merged_fastq: Path to the input merged FASTQ file.
        quality_threshold: The minimum Phred quality score. Reads with 3 consecutive bases below this threshold are discarded.
        min_percent_length: The minimum percentage of the original read length that must remain after trimming.
        qiime2_bin: Path to the QIIME 2 binary (e.g., /opt/qiime2-2023.9/bin).
        filtered_fastq: Path to the output filtered FASTQ file.

    Returns:
        None. The function writes the filtered FASTQ data to the specified output file.

    Raises:
        FileNotFoundError: If the input FASTQ file or QIIME 2 binary is not found.
        ValueError: If the quality_threshold or min_percent_length are invalid.
        subprocess.CalledProcessError: If the QIIME 2 command fails.

    Examples:
        fastq_quality_filter(
            merged_fastq="merged.fastq",
            quality_threshold=20,
            min_percent_length=0.75,
            qiime2_bin="/opt/qiime2-2023.9/bin",
            filtered_fastq="filtered.fastq",
        )
    """

    if not os.path.isfile(merged_fastq):
        raise FileNotFoundError(f"Input FASTQ file not found: {merged_fastq}")
    if not os.path.isdir(qiime2_bin) and not os.path.isfile(qiime2_bin):
        raise FileNotFoundError(f"QIIME 2 binary not found: {qiime2_bin}")

    if not 0 <= min_percent_length <= 1:
        raise ValueError("min_percent_length must be between 0 and 1")

    if quality_threshold < 0:
        raise ValueError("quality_threshold must be non-negative")

    # Construct the QIIME 2 command.
    cmd = [
        os.path.join(qiime2_bin, "qiime"),
        "demux",
        "filter-seqs",
        "--i-seqs",
        merged_fastq,
        "--o-filtered-seqs",
        filtered_fastq,
        "--p-min-quality",
        str(quality_threshold),
        "--p-min-length-fraction",
        str(min_percent_length),
        "--p-max-n",
        "0", # Remove reads with any Ns
    ]

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        raise subprocess.CalledProcessError(
            e.returncode,
            e.cmd,
            output=e.output,
            stderr=e.stderr,
        ) from e

import subprocess
import os
from typing import Optional

def denoise_with_deblur(filtered_fastq: str, qiime2_bin: str, deblur_bin: str, output_biom: str) -> None:
    """Denoises a FASTQ file using Deblur to generate an ASV table in BIOM format.

    Args:
        filtered_fastq: Path to the filtered FASTQ file.
        qiime2_bin: Path to the QIIME 2 bin directory (needed for 'qiime tools export').
        deblur_bin: Path to the Deblur bin directory (needed for 'deblur denoise-16S').
        output_biom: Path to save the resulting ASV table in BIOM format.

    Returns:
        None. The function saves the ASV table to the specified output path.

    Raises:
        FileNotFoundError: If any of the input files or directories do not exist.
        subprocess.CalledProcessError: If any of the subprocess calls fail.
        ValueError: If the input file types are incorrect.

    Examples:
        denoise_with_deblur(
            filtered_fastq='filtered.fastq.gz',
            qiime2_bin='/opt/qiime2-2023.9/bin',
            deblur_bin='/opt/deblur/deblur',
            output_biom='asv_table.biom'
        )
    """

    # Preflight checks
    if not os.path.isfile(filtered_fastq):
        raise FileNotFoundError(f"Filtered FASTQ file not found: {filtered_fastq}")
    if not os.path.isdir(qiime2_bin):
        raise FileNotFoundError(f"QIIME 2 bin directory not found: {qiime2_bin}")
    if not os.path.isdir(deblur_bin):
        raise FileNotFoundError(f"Deblur bin directory not found: {deblur_bin}")

    # Construct paths to executables
    deblur_exe = os.path.join(deblur_bin, 'deblur')
    qiime_tools_export_exe = os.path.join(qiime2_bin, 'qiime')

    if not os.path.isfile(deblur_exe):
        raise FileNotFoundError(f"Deblur executable not found: {deblur_exe}")

    if not os.path.isfile(qiime_tools_export_exe):
        raise FileNotFoundError(f"QIIME tools executable not found: {qiime_tools_export_exe}")

    # Create a temporary directory for intermediate files
    temp_dir = 'temp_deblur_output'
    os.makedirs(temp_dir, exist_ok=True)

    # Deblur denoising
    try:
        subprocess.run(
            [deblur_exe, 'denoise-16S',
             '--seqs-fp', filtered_fastq,
             '--output-dir', temp_dir,
             '--pos-ref-db', os.path.join(deblur_bin, 'core_alignment.fasta'),
             '--neg-ref-db', os.path.join(deblur_bin, 'negative_control_sequences.fasta'),
             '--threads', '1'],  # Consider making threads configurable
            check=True, capture_output=True, text=True
        )
    except subprocess.CalledProcessError as e:
        raise subprocess.CalledProcessError(e.returncode, e.cmd, output=e.output, stderr=e.stderr)

    # Convert to biom format
    try:
        subprocess.run(
            [qiime_tools_export_exe, 'tools', 'export',
             '--input-path', os.path.join(temp_dir, 'deblur.seqs.fa'),
             '--output-path', os.path.join(temp_dir, 'exported_seqs')],
            check=True, capture_output=True, text=True
        )

        # Create the BIOM table using vsearch (or other method)
        with open(os.path.join(temp_dir, 'exported_seqs', 'dna-sequences.fasta')) as fasta_file:
            seqs = fasta_file.readlines()

        # Create feature table
        feature_table = {}
        for i in range(0, len(seqs), 2):
            seq_id = seqs[i].strip()[1:]
            sequence = seqs[i+1].strip()
            feature_table[sequence] = seq_id

        # Write feature table to biom format
        with open(os.path.join(temp_dir, 'feature_table.tsv'), 'w') as outfile:
            outfile.write("#OTU ID\tsequence\n")
            for sequence, seq_id in feature_table.items():
                outfile.write(f"{seq_id}\t{sequence}\n")

        # Convert to biom format using qiime tools import
        subprocess.run(
            [qiime_tools_export_exe, 'tools', 'import',
             '--input-path', os.path.join(temp_dir, 'feature_table.tsv'),
             '--output-path', output_biom,
             '--type', 'FeatureData[Sequence]'],
            check=True, capture_output=True, text=True
        )

    except subprocess.CalledProcessError as e:
        raise subprocess.CalledProcessError(e.returncode, e.cmd, output=e.output, stderr=e.stderr)

    finally:
        # Clean up temporary directory
        import shutil
        shutil.rmtree(temp_dir, ignore_errors=True)

def filter_biom_table(biom_table_path: str, min_reads: int = 1000, output_path: str = "filtered_asv_table.biom") -> None:
    """Filters a BIOM table, removing samples with fewer than a specified number of reads.

    Args:
        biom_table_path: Path to the input BIOM table.
        min_reads: Minimum number of reads for a sample to be retained.
        output_path: Path to save the filtered BIOM table.

    Returns:
        None. The filtered BIOM table is saved to the specified output path.

    Raises:
        FileNotFoundError: If the input BIOM table does not exist.
        ValueError: If the minimum number of reads is not a positive integer.
        ImportError: If the biom-format package is not installed.
        Exception: If there is an error during BIOM table processing.

    Examples:
        filter_biom_table("asv_table.biom", min_reads=1000, output_path="filtered_table.biom")
    """
    try:
        import biom
    except ImportError:
        raise ImportError("The 'biom-format' package is required. Install it with 'pip install biom-format'.")

    if not isinstance(min_reads, int) or min_reads <= 0:
        raise ValueError("Minimum reads must be a positive integer.")

    if not isinstance(biom_table_path, str):
        raise TypeError("biom_table_path must be a string.")

    if not isinstance(output_path, str):
        raise TypeError("output_path must be a string.")

    try:
        with open(biom_table_path) as f:
            pass
    except FileNotFoundError:
        raise FileNotFoundError(f"BIOM table not found at {biom_table_path}")

    try:
        table = biom.load_table(biom_table_path)
        sample_counts = table.sum(axis='sample')
        samples_to_keep = table.ids('sample')[sample_counts >= min_reads]
        filtered_table = table.filter(samples_to_keep, axis='sample')
        with biom.util.biom_open(output_path, 'w') as f:
            filtered_table.to_hdf5(f, "Filtered BIOM table")

    except Exception as e:
        raise Exception(f"Error processing BIOM table: {e}")

import subprocess
import os
from typing import Optional, List

def perform_diversity_analysis(
    asv_table_path: str,
    qiime2_bin: str,
    output_dir: str,
    rarefaction_depth: int = 1000,
    alpha_metrics: Optional[List[str]] = None,
    beta_metric: str = "unweighted_unifrac",
) -> None:
    """Performs diversity analysis on an ASV table using QIIME 2.

    This function executes rarefaction, SRS scaling, and calculates alpha and beta diversity metrics using QIIME 2.

    Args:
        asv_table_path: Path to the filtered ASV table in BIOM format.
        qiime2_bin: Path to the QIIME 2 binary (e.g., /opt/qiime2-2023.9/bin/qiime).
        output_dir: Path to the output directory where diversity results will be stored.
        rarefaction_depth: Rarefaction depth for subsampling.
        alpha_metrics: List of alpha diversity metrics to calculate. If None, defaults to ['observed_features', 'shannon'].
        beta_metric: Beta diversity metric to calculate (e.g., 'unweighted_unifrac', 'weighted_unifrac').

    Returns:
        None. The function writes diversity analysis results to the specified output directory.

    Raises:
        FileNotFoundError: If the ASV table or QIIME 2 binary is not found.
        subprocess.CalledProcessError: If any of the QIIME 2 commands fail.
        ValueError: If the rarefaction depth is invalid.

    Examples:
        perform_diversity_analysis(
            asv_table_path="filtered_asv_table.biom",
            qiime2_bin="/opt/qiime2-2023.9/bin/qiime",
            output_dir="diversity_results",
            rarefaction_depth=1000,
            alpha_metrics=["observed_features", "shannon"],
            beta_metric="unweighted_unifrac",
        )
    """

    if not os.path.exists(asv_table_path):
        raise FileNotFoundError(f"ASV table not found: {asv_table_path}")
    if not os.path.exists(qiime2_bin):
        raise FileNotFoundError(f"QIIME 2 binary not found: {qiime2_bin}")
    if rarefaction_depth <= 0:
        raise ValueError("Rarefaction depth must be a positive integer.")

    os.makedirs(output_dir, exist_ok=True)

    if alpha_metrics is None:
        alpha_metrics = ["observed_features", "shannon"]

    # Rarefaction and SRS scaling
    rarefied_table_path = os.path.join(output_dir, "rarefied_table.biom")
    try:
        subprocess.run(
            [
                qiime2_bin,
                "feature-table",
                "rarefy",
                "--i-table",
                asv_table_path,
                "--p-sampling-depth",
                str(rarefaction_depth),
                "--o-rarefied-table",
                rarefied_table_path,
            ],
            check=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError as e:
        raise subprocess.CalledProcessError(e.returncode, e.cmd, output=e.output, stderr=e.stderr)

    # Alpha diversity
    alpha_vector_path = os.path.join(output_dir, "alpha_vector.qza")
    try:
        subprocess.run(
            [
                qiime2_bin,
                "diversity",
                "alpha",
                "--i-table",
                rarefied_table_path,
                "--p-metric",
                ",".join(alpha_metrics),
                "--o-alpha-diversity",
                alpha_vector_path,
            ],
            check=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError as e:
        raise subprocess.CalledProcessError(e.returncode, e.cmd, output=e.output, stderr=e.stderr)

    # Beta diversity
    beta_diversity_path = os.path.join(output_dir, "beta_diversity.qza")
    try:
        subprocess.run(
            [
                qiime2_bin,
                "diversity",
                "beta",
                "--i-table",
                rarefied_table_path,
                "--p-metric",
                beta_metric,
                "--o-distance-matrix",
                beta_diversity_path,
            ],
            check=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError as e:
        raise subprocess.CalledProcessError(e.returncode, e.cmd, output=e.output, stderr=e.stderr)

    # Visualize alpha diversity
    alpha_visualization_path = os.path.join(output_dir, "alpha_visualization.qzv")
    try:
        subprocess.run(
            [
                qiime2_bin,
                "diversity",
                "alpha-group-significance",
                "--i-alpha-diversity",
                alpha_vector_path,
                "--m-metadata-file",
                asv_table_path, # Using ASV table as metadata (may need a separate metadata file)
                "--o-visualization",
                alpha_visualization_path,
            ],
            check=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError as e:
        print(f"Warning: Alpha group significance visualization failed. Continuing. Error: {e.stderr}")

    # Visualize beta diversity
    beta_visualization_path = os.path.join(output_dir, "beta_visualization.qzv")
    try:
        subprocess.run(
            [
                qiime2_bin,
                "diversity",
                "beta-group-significance",
                "--i-distance-matrix",
                beta_diversity_path,
                "--m-metadata-file",
                asv_table_path, # Using ASV table as metadata (may need a separate metadata file)
                "--o-visualization",
                beta_visualization_path,
                "--p-pairwise",
            ],
            check=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError as e:
        print(f"Warning: Beta group significance visualization failed. Continuing. Error: {e.stderr}")

import subprocess
import os
from typing import List, Tuple


def construct_phylogenetic_tree(
    asv_sequences: List[str],
    mafft_bin: str,
    fasttree_bin: str,
    output_dir: str
) -> Tuple[str, str]:
    """Constructs a phylogenetic tree from ASV sequences using MAFFT and FastTree.

    Args:
        asv_sequences: A list of ASV sequences (strings).  Each sequence should have a unique identifier as its FASTA header.
        mafft_bin: Path to the MAFFT executable.
        fasttree_bin: Path to the FastTree executable.
        output_dir: Directory to store the alignment and tree files.

    Returns:
        A tuple containing the paths to the alignment file (fasta) and the phylogenetic tree file (newick).

    Raises:
        FileNotFoundError: If MAFFT or FastTree executables are not found.
        subprocess.CalledProcessError: If MAFFT or FastTree commands fail.
        ValueError: If the input ASV sequences are invalid.

    Examples:
        >>> asv_seqs = ['>ASV1\nATGC', '>ASV2\nTGCA']
        >>> alignment_file, tree_file = construct_phylogenetic_tree(asv_seqs, '/usr/local/bin/mafft', '/usr/local/bin/FastTree', 'output')
        >>> print(alignment_file)
        output/asv_alignment.fasta
        >>> print(tree_file)
        output/phylo_tree.nwk
    """

    # Preflight checks
    if not os.path.isfile(mafft_bin):
        raise FileNotFoundError(f"MAFFT executable not found at {mafft_bin}")
    if not os.path.isfile(fasttree_bin):
        raise FileNotFoundError(f"FastTree executable not found at {fasttree_bin}")

    if not asv_sequences:
        raise ValueError("Input ASV sequences cannot be empty.")

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    alignment_file = os.path.join(output_dir, "asv_alignment.fasta")
    tree_file = os.path.join(output_dir, "phylo_tree.nwk")

    # Write ASV sequences to a temporary file
    asv_fasta_file = os.path.join(output_dir, "asv_sequences.fasta")
    with open(asv_fasta_file, "w") as f:
        f.write("\n".join(asv_sequences) + "\n")

    # Align sequences using MAFFT
    try:
        subprocess.run(
            [mafft_bin, "--auto", asv_fasta_file],
            check=True,
            capture_output=True,
            text=True,
            stdout=open(alignment_file, "w"),
            stderr=subprocess.PIPE
        )
    except subprocess.CalledProcessError as e:
        raise subprocess.CalledProcessError(e.returncode, e.cmd, output=e.output, stderr=e.stderr) from e

    # Build tree using FastTree
    try:
        subprocess.run(
            [fasttree_bin, "-nt", alignment_file],
            check=True,
            capture_output=True,
            text=True,
            stdout=open(tree_file, "w"),
            stderr=subprocess.PIPE
        )
    except subprocess.CalledProcessError as e:
        raise subprocess.CalledProcessError(e.returncode, e.cmd, output=e.output, stderr=e.stderr) from e

    return alignment_file, tree_file

import subprocess
import os
from typing import Optional

def annotate_asv_with_qiime2(
    asv_fasta: str,
    qiime2_bin: str,
    ezbio_db: str,
    output_tsv: str,
    identity_threshold: float = 0.97
) -> None:
    """Annotates ASV sequences using QIIME2's classify-consensus-vsearch method against the EzBioCloud database.

    Args:
        asv_fasta: Path to the ASV sequences in FASTA format.
        qiime2_bin: Path to the QIIME2 binary (e.g., /opt/qiime2-2023.9/bin/qiime).
        ezbio_db: Path to the EzBioCloud reference database (e.g., a directory containing the .qza).
        output_tsv: Path to the output taxonomy TSV file.
        identity_threshold: Minimum identity threshold for classification (default: 0.97).

    Returns:
        None. Writes the taxonomy to the specified output_tsv file.

    Raises:
        FileNotFoundError: If the ASV FASTA file, QIIME2 binary, or EzBioCloud database are not found.
        ValueError: If the identity threshold is not between 0 and 1.
        subprocess.CalledProcessError: If the QIIME2 command fails.

    Examples:
        annotate_asv_with_qiime2(
            asv_fasta="asv_sequences.fasta",
            qiime2_bin="/opt/qiime2-2023.9/bin/qiime",
            ezbio_db="/path/to/ezbio_db",
            output_tsv="taxonomy.tsv",
            identity_threshold=0.97
        )
    """

    if not os.path.isfile(asv_fasta):
        raise FileNotFoundError(f"ASV FASTA file not found: {asv_fasta}")
    if not os.path.isfile(qiime2_bin):
        raise FileNotFoundError(f"QIIME2 binary not found: {qiime2_bin}")
    if not os.path.exists(ezbio_db):
        raise FileNotFoundError(f"EzBioCloud database not found: {ezbio_db}")

    if not 0 <= identity_threshold <= 1:
        raise ValueError("Identity threshold must be between 0 and 1.")

    cmd = [
        qiime2_bin,
        "feature-classifier",
        "classify-consensus-vsearch",
        "--query-sequences", asv_fasta,
        "--reference-taxonomy", os.path.join(ezbio_db, "taxonomy.qza") if os.path.isdir(ezbio_db) else ezbio_db + '/taxonomy.qza' if os.path.isfile(ezbio_db + '/taxonomy.qza') else ezbio_db,
        "--reference-sequences", os.path.join(ezbio_db, "sequences.qza") if os.path.isdir(ezbio_db) else ezbio_db + '/sequences.qza' if os.path.isfile(ezbio_db + '/sequences.qza') else ezbio_db,
        "--output-taxonomy", output_tsv,
        "--min-id", str(identity_threshold),
        "--threads", "1", # Ensure deterministic behavior
        "--verbose"
    ]

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        raise subprocess.CalledProcessError(e.returncode, e.cmd, e.stdout, e.stderr)
