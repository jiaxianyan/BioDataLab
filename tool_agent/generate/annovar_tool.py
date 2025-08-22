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
