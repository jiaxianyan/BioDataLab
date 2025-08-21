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
