import json
import os
import pickle
import time
from typing import Any

import requests
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from langchain_core.messages import HumanMessage, SystemMessage

import urllib
import subprocess
import shutil
from pathlib import Path
import json
import os
from pathlib import Path
from typing import List, Tuple, Dict, Any, Optional, Literal
import pandas as pd
import GEOparse
from datetime import datetime
from Bio import Entrez, Medline
import time

import openai
from langchain_anthropic import ChatAnthropic

# from langchain_aws import ChatBedrock
from langchain_core.language_models.chat_models import BaseChatModel
from langchain_google_genai import ChatGoogleGenerativeAI
from langchain_ollama import ChatOllama
from langchain_openai import AzureChatOpenAI, ChatOpenAI

import re


def pubmed_fetch_abstracts(
    query: str,
    date_from: str,                 # e.g. "2020/01/01"
    date_to: str,                   # e.g. "2025/12/31"
    email: str,                     # contact email required by NCBI
    out_dir: str = "benchmark/results/3",
    out_tsv: Optional[str] = None,  # default: <out_dir>/pubmed_abstracts.tsv
    api_key: Optional[str] = None,
    retmax: int = 10000,            # esearch page size
    batch_size: int = 200,          # efetch batch size (<=200 recommended)
    sleep_between: float = 0.34     # polite delay between requests
) -> Tuple[List[str], Path]:
    """
    Search PubMed with a keyword query and publication date range, retrieve PMIDs via
    esearch, then fetch titles and abstracts via efetch (MEDLINE). Writes a TSV with
    columns: pmid, title, abstract.

    Parameters
    ----------
    query : str
        PubMed E-utilities query string, e.g. '("SARS-CoV-2" AND "immune escape")'.
    date_from : str
        Lower bound for publication date filter (YYYY/MM/DD), applied to PDAT.
    date_to : str
        Upper bound for publication date filter (YYYY/MM/DD), applied to PDAT.
    email : str
        Contact email required by NCBI E-utilities.
    out_dir : str
        Output directory to store the TSV file.
    out_tsv : Optional[str]
        Explicit output TSV path; if None, defaults to '<out_dir>/pubmed_abstracts.tsv'.
    api_key : Optional[str]
        NCBI API key to increase rate limits (recommended).
    retmax : int
        Page size for esearch pagination.
    batch_size : int
        Batch size for efetch requests (<=200 is recommended).
    sleep_between : float
        Delay seconds between NCBI requests to avoid rate limiting.

    Returns
    -------
    pmid_list : List[str]
        List of PMIDs retrieved from PubMed (unique, order-preserving).
    abstracts_path : Path
        Path to the TSV with columns [pmid, title, abstract].

    Notes
    -----
    - Uses PDAT as the publication date field. Adjust as needed.
    - Results can drift over time as PubMed updates; freeze the TSV as refdata for benchmarks.
    """
    # --- Init Entrez ---
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    # --- Compose query with date range on PDAT ---
    df = date_from or "1800/01/01"
    dt = date_to or datetime.today().strftime("%Y/%m/%d")
    term = f'({query}) AND ("{df}"[PDAT] : "{dt}"[PDAT])'

    out_root = Path(out_dir)
    out_root.mkdir(parents=True, exist_ok=True)
    abstracts_path = Path(out_tsv) if out_tsv else (out_root / "pubmed_abstracts.tsv")

    # --- 1) esearch (paginated) to collect all PMIDs ---
    pmids: List[str] = []
    retstart = 0
    while True:
        with Entrez.esearch(db="pubmed", term=term, retmax=retmax, retstart=retstart) as h:
            rec = Entrez.read(h)
        ids = rec.get("IdList", [])
        pmids.extend(ids)
        retstart += len(ids)
        time.sleep(sleep_between)
        if len(ids) < retmax:
            break

    # Deduplicate while preserving order
    seen = set()
    pmids = [p for p in pmids if not (p in seen or seen.add(p))]

    # Early exit if empty
    if not pmids:
        # Write empty TSV with header for consistency
        pd.DataFrame(columns=["pmid", "title", "abstract"]).to_csv(abstracts_path, sep="\t", index=False)
        return [], abstracts_path

    # --- 2) efetch in batches (MEDLINE) ---
    rows = []
    for i in range(0, len(pmids), batch_size):
        chunk = pmids[i:i+batch_size]
        with Entrez.efetch(db="pubmed", id=",".join(chunk), rettype="medline", retmode="text") as h:
            records = list(Medline.parse(h))
        for r in records:
            pmid = str(r.get("PMID", "")).strip()
            title = (r.get("TI") or "").strip()
            abstract = (r.get("AB") or "").strip()
            if pmid:
                rows.append({"pmid": pmid, "title": title, "abstract": abstract})
        time.sleep(sleep_between)

    # --- 3) Write TSV (sorted by PMID for stability) ---
    df_out = pd.DataFrame(rows).drop_duplicates(subset=["pmid"]).sort_values("pmid")
    df_out.to_csv(abstracts_path, sep="\t", index=False)

    return [str(p) for p in df_out["pmid"].tolist()], abstracts_path



def download_sra_runs(query_result: dict, sra_toolkit_dir: str = "./operation_env/tool_lake/sra_toolkit", outdir: str = "./sra_fastq"):
    """
    Download SRA runs as FASTQ files using SRA Toolkit, based on query_sra results.

    Parameters
    ----------
    query_result : dict
        Result dictionary from query_sra(). Must contain "formatted_results" with SRA Run IDs.
    sra_toolkit_dir : str
        Path to SRA Toolkit binaries (must contain prefetch/fasterq-dump).
        Example: "/home/user/sratoolkit/bin"
    outdir : str, default="fastq"
        Directory to save downloaded FASTQ files.

    Returns
    -------
    list
        List of paths to downloaded FASTQ files.

    Raises
    ------
    RuntimeError
        If query_result is invalid or SRA Toolkit commands fail.

    Example
    -------
    result = query_sra(prompt="breast cancer RNA-seq")
    download_sra_runs(result, "/usr/local/sratoolkit/bin", outdir="data")
    """
    if not query_result or "formatted_results" not in query_result:
        raise RuntimeError("Invalid query_result: missing formatted_results")

    # Ensure output directory exists
    os.makedirs(outdir, exist_ok=True)

    # Extract run IDs (SRRxxxxxx)
    run_ids = []
    for entry in query_result["formatted_results"]:
        if isinstance(entry, dict):
            # Common field in NCBI esummary output for SRA
            if "Runs" in entry:
                run_ids.extend([r["acc"] for r in entry["Runs"] if "acc" in r])
            elif "uid" in entry and entry["uid"].startswith("SRR"):
                run_ids.append(entry["uid"])
        elif isinstance(entry, str) and entry.startswith("SRR"):
            run_ids.append(entry)

    if not run_ids:
        raise RuntimeError("No SRR run IDs found in query_result")

    # Paths to SRA Toolkit binaries
    prefetch = os.path.join(sra_toolkit_dir, "prefetch")
    fasterq_dump = os.path.join(sra_toolkit_dir, "fasterq-dump")

    downloaded_files = []

    for run_id in run_ids:
        print(f"[INFO] Downloading {run_id}...")
        try:
            # Step 1: Prefetch SRA file
            subprocess.run([prefetch, run_id], check=True)

            # Step 2: Convert to FASTQ
            subprocess.run(
                [fasterq_dump, run_id, "-O", outdir, "--split-files"],
                check=True
            )

            # Collect generated FASTQ files
            fq1 = os.path.join(outdir, f"{run_id}_1.fastq")
            fq2 = os.path.join(outdir, f"{run_id}_2.fastq")
            if os.path.exists(fq1):
                downloaded_files.append(fq1)
            if os.path.exists(fq2):
                downloaded_files.append(fq2)

        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"SRA Toolkit failed for {run_id}: {str(e)}")

    return downloaded_files

def download_geo(query_result: dict, outdir: str = "./geo_data", file_type: str = "suppl"):
    """
    Download GEO data (GSE/GSM) from NCBI FTP based on query_geo results.

    Parameters
    ----------
    query_result : dict
        Result dictionary from query_geo(). Must contain "formatted_results" with GEO IDs (GSE, GSM).
    outdir : str, default="geo_data"
        Directory to save downloaded files.
    file_type : str, default="suppl"
        Which subdirectory to download from (suppl, soft, miniml).

    Returns
    -------
    list
        List of downloaded file paths.

    Raises
    ------
    RuntimeError
        If no GEO IDs found or download fails.

    Example
    -------
    result = query_geo(prompt="breast cancer RNA-seq")
    files = download_geo(result, outdir="data", file_type="suppl")
    """
    if not query_result or "formatted_results" not in query_result:
        raise RuntimeError("Invalid query_result: missing formatted_results")

    os.makedirs(outdir, exist_ok=True)

    # Extract GEO IDs (GSE / GSM)
    geo_ids = []
    for entry in query_result["formatted_results"]:
        if isinstance(entry, dict):
            if "uid" in entry and re.match(r"^GS[EM]\d+$", entry["uid"]):
                geo_ids.append(entry["uid"])
            elif "Accession" in entry and re.match(r"^GS[EM]\d+$", entry["Accession"]):
                geo_ids.append(entry["Accession"])
        elif isinstance(entry, str) and re.match(r"^GS[EM]\d+$", entry):
            geo_ids.append(entry)

    if not geo_ids:
        raise RuntimeError("No GEO IDs (GSE/GSM) found in query_result")

    downloaded_files = []

    for geo_id in geo_ids:
        # Build FTP URL (grouping rule: GSE12345 → GSE12nnn/GSE12345/)
        prefix = geo_id[:-3] + "nnn"  # e.g. GSE12345 → GSE12nnn
        url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{geo_id}/{file_type}/"

        print(f"[INFO] Checking GEO FTP: {url}")

        try:
            with urllib.request.urlopen(url) as resp:
                html = resp.read().decode("utf-8")
        except Exception:
            print(f"[WARN] Cannot access {url}")
            continue

        # Extract file links
        files = re.findall(r'href="([^"]+)"', html)
        files = [f for f in files if not f.startswith("?") and not f.endswith("/")]

        for fname in files:
            f_url = url + fname
            f_out = os.path.join(outdir, fname)
            try:
                print(f"[INFO] Downloading {f_url} → {f_out}")
                urllib.request.urlretrieve(f_url, f_out)
                downloaded_files.append(f_out)
            except Exception as e:
                print(f"[ERROR] Failed to download {f_url}: {e}")

    return downloaded_files