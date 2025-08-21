

import subprocess
import os
from typing import Optional, List, Dict


def run_trimmomatic(
    trimmomatic_path: str,
    mode: str,
    input_r1: Optional[str] = None,
    input_r2: Optional[str] = None,
    output_r1_paired: Optional[str] = None,
    output_r1_unpaired: Optional[str] = None,
    output_r2_paired: Optional[str] = None,
    output_r2_unpaired: Optional[str] = None,
    input_fastq: Optional[str] = None,
    output_fastq: Optional[str] = None,
    adapters_path: Optional[str] = None,
    leading: Optional[int] = None,
    trailing: Optional[int] = None,
    slidingwindow: Optional[str] = None,
    minlen: Optional[int] = None,
    extra_args: Optional[List[str]] = None,
) -> Dict:
    """Runs Trimmomatic for quality trimming and adapter removal.

    Args:
        trimmomatic_path: The full path to the Trimmomatic executable (e.g., /path/to/trimmomatic.jar).
        mode: 'PE' for paired-end or 'SE' for single-end.
        input_r1: Path to the input R1 FASTQ file (required for PE).
        input_r2: Path to the input R2 FASTQ file (required for PE).
        output_r1_paired: Path to the output R1 paired FASTQ file (required for PE).
        output_r1_unpaired: Path to the output R1 unpaired FASTQ file (required for PE).
        output_r2_paired: Path to the output R2 paired FASTQ file (required for PE).
        output_r2_unpaired: Path to the output R2 unpaired FASTQ file (required for PE).
        input_fastq: Path to the input FASTQ file (required for SE).
        output_fastq: Path to the output FASTQ file (required for SE).
        adapters_path: Path to the adapters FASTA file.
        leading: Minimum quality required to keep a leading base.
        trailing: Minimum quality required to keep a trailing base.
        slidingwindow: Sliding window size and minimum quality (e.g., "4:20").
        minlen: Minimum read length to keep.
        extra_args: A list of extra arguments to pass to Trimmomatic.

    Returns:
        A dictionary containing the return code, stdout, stderr, and output file paths.
        Example: {"returncode": 0, "stdout": "...", "stderr": "...", "outputs": {"output_fastq": "/path/to/output.fastq"}}

    Raises:
        FileNotFoundError: If the Trimmomatic executable or input files are not found.
        ValueError: If the mode is invalid or required input/output files are missing.
        OSError: If the output directory is not writable.

    Examples:
        # Paired-end mode
        result = run_trimmomatic(
            trimmomatic_path="/path/to/trimmomatic.jar",
            mode="PE",
            input_r1="/path/to/input_R1.fastq",
            input_r2="/path/to/input_R2.fastq",
            output_r1_paired="/path/to/output_R1_paired.fastq",
            output_r1_unpaired="/path/to/output_R1_unpaired.fastq",
            output_r2_paired="/path/to/output_R2_paired.fastq",
            output_r2_unpaired="/path/to/output_R2_unpaired.fastq",
            adapters_path="/path/to/adapters.fasta",
            leading=3,
            trailing=3,
            slidingwindow="4:15",
            minlen=36,
        )

        # Single-end mode
        result = run_trimmomatic(
            trimmomatic_path="/path/to/trimmomatic.jar",
            mode="SE",
            input_fastq="/path/to/input.fastq",
            output_fastq="/path/to/output.fastq",
            leading=3,
            trailing=3,
            slidingwindow="4:15",
            minlen=36,
        )
    """

    # Preflight validations
    if not os.path.exists(trimmomatic_path):
        raise FileNotFoundError(f"Trimmomatic executable not found: {trimmomatic_path}")
    if not os.access(trimmomatic_path, os.X_OK) and not trimmomatic_path.endswith(".jar"):
        raise OSError(f"Trimmomatic executable not executable: {trimmomatic_path}")

    command = []
    if trimmomatic_path.endswith(".jar"):
      command.extend(['java', '-jar', trimmomatic_path])
    else:
      command.append(trimmomatic_path)

    if mode not in ("PE", "SE"):
        raise ValueError(f"Invalid mode: {mode}. Must be 'PE' or 'SE'.")
    command.append(mode)

    if mode == "PE":
        if not all([input_r1, input_r2, output_r1_paired, output_r1_unpaired, output_r2_paired, output_r2_unpaired]):
            raise ValueError("Missing required input/output files for PE mode.")
        if not os.path.exists(input_r1):
            raise FileNotFoundError(f"Input R1 file not found: {input_r1}")
        if not os.path.exists(input_r2):
            raise FileNotFoundError(f"Input R2 file not found: {input_r2}")

        command.extend([input_r1, input_r2, output_r1_paired, output_r1_unpaired, output_r2_paired, output_r2_unpaired])
        outputs = {
            "output_r1_paired": output_r1_paired,
            "output_r1_unpaired": output_r1_unpaired,
            "output_r2_paired": output_r2_paired,
            "output_r2_unpaired": output_r2_unpaired,
        }

    elif mode == "SE":
        if not all([input_fastq, output_fastq]):
            raise ValueError("Missing required input/output files for SE mode.")
        if not os.path.exists(input_fastq):
            raise FileNotFoundError(f"Input FASTQ file not found: {input_fastq}")
        outputs = {"output_fastq": output_fastq}
        command.extend([input_fastq, output_fastq])

    if adapters_path:
        if not os.path.exists(adapters_path):
            raise FileNotFoundError(f"Adapters file not found: {adapters_path}")
        command.append(f"ILLUMINACLIP:{adapters_path}:2:30:10")

    if leading is not None:
        command.append(f"LEADING:{leading}")
    if trailing is not None:
        command.append(f"TRAILING:{trailing}")
    if slidingwindow is not None:
        command.append(f"SLIDINGWINDOW:{slidingwindow}")
    if minlen is not None:
        command.append(f"MINLEN:{minlen}")

    if extra_args:
        command.extend(extra_args)

    try:
        result = subprocess.run(command, capture_output=True, text=True, check=False)
        return {
            "returncode": result.returncode,
            "stdout": result.stdout,
            "stderr": result.stderr,
            "outputs": outputs,
        }
    except Exception as e:
        return {
            "returncode": 1,
            "stdout": "",
            "stderr": str(e),
            "outputs": {},
        }

import subprocess
import os
from typing import List, Optional, Dict


def run_fastqc(
    fastqc_executable: str,
    input_files: List[str],
    output_dir: str,
    threads: int = 1,
    extra_args: Optional[List[str]] = None
) -> Dict:
    """Runs FastQC on the given input files and returns the execution results.

    Args:
        fastqc_executable: The full path to the FastQC executable.
        input_files: A list of input FASTQ or BAM files.
        output_dir: The directory where FastQC output will be written.
        threads: The number of threads to use.
        extra_args: A list of additional arguments to pass to FastQC.

    Returns:
        A dictionary containing the return code, stdout, stderr, and output file paths.

    Raises:
        FileNotFoundError: If the FastQC executable or any input file does not exist.
        OSError: If the output directory cannot be created or is not writable.
        ValueError: If the number of threads is not a positive integer.

    Examples:
        >>> result = run_fastqc(
        ...     fastqc_executable='/path/to/fastqc',
        ...     input_files=['/path/to/input1.fastq', '/path/to/input2.fastq'],
        ...     output_dir='/path/to/output',
        ...     threads=4
        ... )
        >>> print(result['returncode'])
        0
        >>> print(result['stdout'])
        'FastQC completed successfully...'
    """

    # Preflight validations
    if not os.path.exists(fastqc_executable) or not os.access(fastqc_executable, os.X_OK):
        raise FileNotFoundError(f"FastQC executable not found or not executable: {fastqc_executable}")

    for input_file in input_files:
        if not os.path.exists(input_file):
            raise FileNotFoundError(f"Input file not found: {input_file}")

    if not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir, exist_ok=True)
        except OSError as e:
            raise OSError(f"Cannot create output directory: {output_dir} - {e}")

    if not os.access(output_dir, os.W_OK):
        raise OSError(f"Output directory not writable: {output_dir}")

    if not isinstance(threads, int) or threads <= 0:
        raise ValueError("Number of threads must be a positive integer.")

    # Build command
    command = [
        fastqc_executable,
        "--threads",
        str(threads),
        "--outdir",
        output_dir,
        *input_files,
    ]

    if extra_args:
        command.extend(extra_args)

    # Execute command
    try:
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=False  # Do not raise exceptions on non-zero return codes
        )

        # Determine output files (assuming FastQC names them predictably)
        outputs = {}
        for input_file in input_files:
            base_name = os.path.basename(input_file)
            outputs[input_file] = os.path.join(output_dir, f"{os.path.splitext(base_name)[0]}_fastqc.html")

        return {
            "returncode": result.returncode,
            "stdout": result.stdout,
            "stderr": result.stderr,
            "outputs": outputs
        }
    except Exception as e:
        return {
            "returncode": 1,
            "stdout": "",
            "stderr": str(e),
            "outputs": {}
        }
