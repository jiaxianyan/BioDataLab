import os
import json
import shlex
import subprocess
import datetime
from typing import Optional, Dict, Any, List, Tuple, Sequence, Union
import tempfile
import shutil


def _ensure_file(path: str, desc: str):
    if not os.path.isfile(path):
        raise FileNotFoundError(f"{desc} not found: {path}")

def _run(cmd: Sequence[str], cwd: Optional[str]=None) -> subprocess.CompletedProcess:
    # 不依赖 PATH：cmd[0] 必须是绝对路径或相对可执行路径
    res = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, check=False)
    if res.returncode != 0:
        raise RuntimeError(
            "Command failed:\n"
            f"$ {' '.join(shlex.quote(c) for c in cmd)}\n"
            f"STDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}"
        )
    return res

def _ensure_dir(p: str) -> str:
    os.makedirs(p, exist_ok=True)
    return os.path.abspath(p)

def _touch_parent(p: str) -> None:
    _ensure_dir(os.path.dirname(p))

def _is_fresh(output_path: str, input_paths: List[str]) -> bool:
    if not os.path.exists(output_path):
        return False
    out_mtime = os.path.getmtime(output_path)
    for ip in input_paths:
        if ip and os.path.exists(ip) and os.path.getmtime(ip) > out_mtime:
            return False
    return True

def _run_qiime(cmd: List[str], cwd: Optional[str] = None, log_file: Optional[str] = None) -> None:
    try:
        proc = subprocess.run(cmd, cwd=cwd, check=True, capture_output=True, text=True)
        if log_file:
            with open(log_file, "a", encoding="utf-8") as f:
                f.write(f"[{datetime.datetime.now().isoformat()}] CMD: {' '.join(shlex.quote(c) for c in cmd)}\n")
                f.write(proc.stdout + "\n")
                f.write(proc.stderr + "\n")
    except subprocess.CalledProcessError as e:
        msg = f"Command failed with exit code {e.returncode}: {' '.join(shlex.quote(c) for c in cmd)}\nSTDOUT:\n{e.stdout}\nSTDERR:\n{e.stderr}"
        if log_file:
            with open(log_file, "a", encoding="utf-8") as f:
                f.write(f"[{datetime.datetime.now().isoformat()}] ERROR: {msg}\n")
        raise RuntimeError(msg) from e

def _require_file(p: Optional[str], name: str) -> str:
    if not p or not os.path.exists(p):
        raise FileNotFoundError(f"Missing required file for '{name}': {p}")
    return os.path.abspath(p)

def _provenance(output_dir: str, info: Dict[str, Any]) -> None:
    path = os.path.join(output_dir, "provenance.jsonl")
    with open(path, "a", encoding="utf-8") as f:
        f.write(json.dumps({"ts": datetime.datetime.now().isoformat(), **info}, ensure_ascii=False) + "\n")

# 1) 数据导入与初步质控
def q2_ingest_and_qc(
    manifest_csv: str,
    output_dir: str,
    qiime_bin: str = "./operation_env/tool_lake/qiime",
    sample_metadata_tsv: Optional[str] = None,
    paired_end: bool = True,
    join_paired: bool = False,
    cutadapt_params: Optional[str] = None,
    qscore_filter: bool = False,
    trunc_len_f: Optional[int] = None,
    trunc_len_r: Optional[int] = None,
    phred_offset: int = 33,
    force: bool = False,
) -> Dict[str, str]:
    """
    Import FASTQ via manifest, optionally join paired reads, trim adapters, perform q-score filtering,
    and summarize demultiplexed data. Returns demux artifact and visualization paths.
    """
    manifest_csv = _require_file(manifest_csv, "manifest_csv")
    if sample_metadata_tsv:
        _require_file(sample_metadata_tsv, "sample_metadata_tsv")
    qiime = _require_file(qiime_bin, "qiime_bin")
    out = _ensure_dir(output_dir)
    log_path = os.path.join(out, "run.log")

    fmt = "PairedEndFastqManifestPhred{}V2".format(phred_offset) if paired_end else "SingleEndFastqManifestPhred{}V2".format(phred_offset)
    demux_qza = os.path.join(out, "demux.qza")
    demux_qzv = os.path.join(out, "demux.qzv")

    if force or not os.path.exists(demux_qza):
        _run_qiime([qiime, "tools", "import",
              "--type", "SampleData[PairedEndSequencesWithQuality]" if paired_end else "SampleData[SequencesWithQuality]",
              "--input-path", manifest_csv,
              "--input-format", fmt,
              "--output-path", demux_qza], log_file=log_path)

    work_qza = demux_qza

    if paired_end and join_paired:
        joined_qza = os.path.join(out, "joined.qza")
        joined_qzv = os.path.join(out, "joined.qzv")
        if force or not os.path.exists(joined_qza):
            _run_qiime([qiime, "vsearch", "join-pairs",
                  "--i-demultiplexed-seqs", demux_qza,
                  "--o-joined-sequences", joined_qza,
                  "--o-join-summary", joined_qzv], log_file=log_path)
        work_qza = joined_qza

    if cutadapt_params:
        trimmed_qza = os.path.join(out, "trimmed.qza")
        # user passes raw parameters string for cutadapt (e.g., "--p-front-f ACTG --p-front-r ACTG")
        if force or not os.path.exists(trimmed_qza):
            cmd = [qiime, "cutadapt", "trim-paired" if paired_end else "trim-single",
                   "--i-demultiplexed-sequences", work_qza,
                   "--o-trimmed-sequences", trimmed_qza]
            cmd.extend(shlex.split(cutadapt_params))
            _run_qiime(cmd, log_file=log_path)
        work_qza = trimmed_qza

    if qscore_filter:
        filt_qza = os.path.join(out, "qscore_filtered.qza")
        filt_stats_qzv = os.path.join(out, "qscore_stats.qzv")
        if force or not os.path.exists(filt_qza):
            _run_qiime([qiime, "quality-filter", "q-score-joined" if (paired_end and join_paired) else "q-score",
                  "--i-demux", work_qza,
                  "--o-filtered-sequences", filt_qza,
                  "--o-filter-stats", filt_stats_qzv], log_file=log_path)
        work_qza = filt_qza

    if paired_end and (trunc_len_f or trunc_len_r):
        # Optional truncation via dada2 pre-truncation is generally applied inside denoise;
        # here we just record intent in provenance to keep interface simple.
        _provenance(out, {"note": "Truncation parameters recorded; actual truncation should be done in denoise step.",
                          "trunc_len_f": trunc_len_f, "trunc_len_r": trunc_len_r})

    if force or not os.path.exists(demux_qzv) or not _is_fresh(demux_qzv, [work_qza]):
        _run_qiime([qiime, "demux", "summarize",
              "--i-data", work_qza,
              "--o-visualization", demux_qzv], log_file=log_path)

    _provenance(out, {"stage": "ingest_and_qc", "manifest": os.path.abspath(manifest_csv), "paired_end": paired_end,
                      "join_paired": join_paired, "qscore_filter": qscore_filter, "cutadapt_params": cutadapt_params})

    return {"demux_qza": os.path.abspath(work_qza), "demux_qzv": os.path.abspath(demux_qzv)}

# 2) 去噪与特征表
def q2_denoise_and_feature_table(
    demux_qza: str,
    output_dir: str,
    qiime_bin: str = "./operation_env/tool_lake/qiime",
    method: str = "deblur",  # 'deblur' or 'dada2'
    trim_length: Optional[int] = 250,  # for deblur
    trunc_len_f: Optional[int] = None,  # for dada2 paired
    trunc_len_r: Optional[int] = None,  # for dada2 paired
    min_reads_per_sample: Optional[int] = None,
    threads: int = 0,
    force: bool = False,
) -> Dict[str, str]:
    """
    Denoise reads (Deblur or DADA2), generate feature table and representative sequences,
    and optionally filter samples by minimal reads.
    """
    demux_qza = _require_file(demux_qza, "demux_qza")
    qiime = _require_file(qiime_bin, "qiime_bin")
    out = _ensure_dir(output_dir)
    log_path = os.path.join(out, "run.log")

    table_qza = os.path.join(out, "table.qza")
    rep_qza = os.path.join(out, "rep-seqs.qza")
    stats_qzv = os.path.join(out, f"{method}_stats.qzv")

    if method.lower() == "deblur":
        if force or not os.path.exists(table_qza) or not os.path.exists(rep_qza):
            cmd = [qiime, "deblur", "denoise-16S",
                   "--i-demultiplexed-seqs", demux_qza,
                   "--o-representative-sequences", rep_qza,
                   "--o-table", table_qza,
                   "--o-stats", stats_qzv]
            if trim_length:
                cmd += ["--p-trim-length", str(trim_length)]
            if threads and threads > 0:
                cmd += ["--p-jobs-to-start", str(threads)]
            _run_qiime(cmd, log_file=log_path)
    elif method.lower() == "dada2":
        # Detect single vs paired via filename convention is unreliable; expect user to choose trunc params accordingly.
        if trunc_len_f is None and trunc_len_r is None:
            # Assume single-end
            if force or not os.path.exists(table_qza) or not os.path.exists(rep_qza):
                cmd = [qiime, "dada2", "denoise-single",
                       "--i-demultiplexed-seqs", demux_qza,
                       "--o-representative-sequences", rep_qza,
                       "--o-table", table_qza,
                       "--o-denoising-stats", stats_qzv]
                if trunc_len_f:
                    cmd += ["--p-trunc-len", str(trunc_len_f)]
                if threads and threads > 0:
                    cmd += ["--p-n-threads", str(threads)]
                _run_qiime(cmd, log_file=log_path)
        else:
            # Paired-end
            if force or not os.path.exists(table_qza) or not os.path.exists(rep_qza):
                cmd = [qiime, "dada2", "denoise-paired",
                       "--i-demultiplexed-seqs", demux_qza,
                       "--o-representative-sequences", rep_qza,
                       "--o-table", table_qza,
                       "--o-denoising-stats", stats_qzv]
                if trunc_len_f:
                    cmd += ["--p-trunc-len-f", str(trunc_len_f)]
                if trunc_len_r:
                    cmd += ["--p-trunc-len-r", str(trunc_len_r)]
                if threads and threads > 0:
                    cmd += ["--p-n-threads", str(threads)]
                _run_qiime(cmd, log_file=log_path)
    else:
        raise ValueError("method must be 'deblur' or 'dada2'.")

    filtered_table_qza = ""
    if min_reads_per_sample and min_reads_per_sample > 0:
        filtered_table_qza = os.path.join(out, f"table_min{min_reads_per_sample}.qza")
        if force or not os.path.exists(filtered_table_qza) or not _is_fresh(filtered_table_qza, [table_qza]):
            _run_qiime([qiime, "feature-table", "filter-samples",
                  "--i-table", table_qza,
                  "--p-min-frequency", str(min_reads_per_sample),
                  "--o-filtered-table", filtered_table_qza], log_file=log_path)

    _provenance(out, {"stage": "denoise_and_feature_table", "method": method,
                      "trim_length": trim_length, "trunc_len_f": trunc_len_f, "trunc_len_r": trunc_len_r,
                      "min_reads_per_sample": min_reads_per_sample})

    result = {
        "table_qza": os.path.abspath(filtered_table_qza or table_qza),
        "rep_seqs_qza": os.path.abspath(rep_qza),
        "denoise_stats_qzv": os.path.abspath(stats_qzv),
    }
    if filtered_table_qza:
        result["table_filtered_qza"] = os.path.abspath(filtered_table_qza)
    return result

# 3) 分类学注释
def q2_taxonomy(
    table_qza: str,
    rep_seqs_qza: str,
    output_dir: str,
    qiime_bin: str = "./operation_env/tool_lake/qiime",
    metadata_tsv: Optional[str] = None,
    pretrained_classifier_qza: Optional[str] = None,
    ref_seqs_qza: Optional[str] = None,
    ref_taxonomy_qza: Optional[str] = None,
    method: str = "sklearn",  # 'sklearn' or 'vsearch'
    perc_identity: float = 0.97,
    force: bool = False,
) -> Dict[str, str]:
    """
    Assign taxonomy using a pretrained classifier (preferred) or by training from reference,
    then generate taxa barplot.
    """
    table_qza = _require_file(table_qza, "table_qza")
    rep_seqs_qza = _require_file(rep_seqs_qza, "rep_seqs_qza")
    qiime = _require_file(qiime_bin, "qiime_bin")
    if metadata_tsv:
        _require_file(metadata_tsv, "metadata_tsv")
    out = _ensure_dir(output_dir)
    log_path = os.path.join(out, "run.log")

    taxonomy_qza = os.path.join(out, "taxonomy.qza")
    taxonomy_qzv = os.path.join(out, "taxonomy.qzv")

    classifier_qza = None
    if pretrained_classifier_qza:
        classifier_qza = _require_file(pretrained_classifier_qza, "pretrained_classifier_qza")
    else:
        # Need to train classifier if using sklearn without pretrained, or prepare for vsearch consensus
        if not (ref_seqs_qza and ref_taxonomy_qza):
            raise ValueError("Either provide 'pretrained_classifier_qza' OR both 'ref_seqs_qza' and 'ref_taxonomy_qza'.")
        ref_seqs_qza = _require_file(ref_seqs_qza, "ref_seqs_qza")
        ref_taxonomy_qza = _require_file(ref_taxonomy_qza, "ref_taxonomy_qza")
        if method == "sklearn":
            classifier_qza = os.path.join(out, "trained_classifier.qza")
            if force or not os.path.exists(classifier_qza):
                _run_qiime([qiime, "feature-classifier", "fit-classifier-naive-bayes",
                      "--i-reference-reads", ref_seqs_qza,
                      "--i-reference-taxonomy", ref_taxonomy_qza,
                      "--o-classifier", classifier_qza], log_file=log_path)

    if method == "sklearn":
        if force or not os.path.exists(taxonomy_qza):
            _run_qiime([qiime, "feature-classifier", "classify-sklearn",
                  "--i-classifier", classifier_qza,
                  "--i-reads", rep_seqs_qza,
                  "--o-classification", taxonomy_qza], log_file=log_path)
    elif method == "vsearch":
        if force or not os.path.exists(taxonomy_qza):
            if not (ref_seqs_qza and ref_taxonomy_qza):
                raise ValueError("VSEARCH method requires 'ref_seqs_qza' and 'ref_taxonomy_qza'.")
            _run_qiime([qiime, "feature-classifier", "classify-consensus-vsearch",
                  "--i-query", rep_seqs_qza,
                  "--i-reference-reads", ref_seqs_qza,
                  "--i-reference-taxonomy", ref_taxonomy_qza,
                  "--p-perc-identity", str(perc_identity),
                  "--o-classification", taxonomy_qza], log_file=log_path)
    else:
        raise ValueError("method must be 'sklearn' or 'vsearch'.")

    taxa_barplot_qzv = os.path.join(out, "taxa-barplot.qzv")
    if force or not os.path.exists(taxa_barplot_qzv):
        cmd = [qiime, "taxa", "barplot",
               "--i-table", table_qza,
               "--i-taxonomy", taxonomy_qza,
               "--o-visualization", taxa_barplot_qzv]
        if metadata_tsv:
            cmd += ["--m-metadata-file", metadata_tsv]
        _run_qiime(cmd, log_file=log_path)

    _provenance(out, {"stage": "taxonomy", "method": method, "perc_identity": perc_identity,
                      "used_pretrained": bool(pretrained_classifier_qza)})

    result = {
        "taxonomy_qza": os.path.abspath(taxonomy_qza),
        "taxonomy_qzv": os.path.abspath(taxonomy_qzv),
        "taxa_barplot_qzv": os.path.abspath(taxa_barplot_qzv),
    }
    if pretrained_classifier_qza is None and method == "sklearn":
        result["trained_classifier_qza"] = os.path.abspath(classifier_qza)
    return result

# 4) 系统发育与核心多样性
def q2_phylogeny_and_core_diversity(
    table_qza: str,
    rep_seqs_qza: str,
    output_dir: str,
    qiime_bin: str = "./operation_env/tool_lake/qiime",
    sampling_depth: int = 1000,
    use_midpoint_root: bool = True,
    force: bool = False,
) -> Dict[str, str]:
    """
    Build alignment and phylogenetic tree, then run core-metrics-phylogenetic to compute
    alpha/beta diversity, PCoA, and Emperor visualizations.
    """
    table_qza = _require_file(table_qza, "table_qza")
    rep_seqs_qza = _require_file(rep_seqs_qza, "rep_seqs_qza")
    qiime = _require_file(qiime_bin, "qiime_bin")
    out = _ensure_dir(output_dir)
    log_path = os.path.join(out, "run.log")

    aligned_qza = os.path.join(out, "aligned-rep-seqs.qza")
    masked_qza = os.path.join(out, "masked-alignment.qza")
    tree_qza = os.path.join(out, "unrooted-tree.qza")
    rooted_tree_qza = os.path.join(out, "rooted-tree.qza")

    if force or not os.path.exists(aligned_qza):
        _run_qiime([qiime, "alignment", "mafft",
              "--i-sequences", rep_seqs_qza,
              "--o-alignment", aligned_qza], log_file=log_path)

    if force or not os.path.exists(masked_qza):
        _run_qiime([qiime, "alignment", "mask",
              "--i-alignment", aligned_qza,
              "--o-masked-alignment", masked_qza], log_file=log_path)

    if force or not os.path.exists(tree_qza):
        _run_qiime([qiime, "phylogeny", "fasttree",
              "--i-alignment", masked_qza,
              "--o-tree", tree_qza], log_file=log_path)

    if force or not os.path.exists(rooted_tree_qza):
        if use_midpoint_root:
            _run_qiime([qiime, "phylogeny", "midpoint-root",
                  "--i-tree", tree_qza,
                  "--o-rooted-tree", rooted_tree_qza], log_file=log_path)
        else:
            _run_qiime([qiime, "phylogeny", "root",
                  "--i-tree", tree_qza,
                  "--o-rooted-tree", rooted_tree_qza], log_file=log_path)

    core_dir = os.path.join(out, "core-metrics")
    _ensure_dir(core_dir)
    # core-metrics-phylogenetic will create a directory with multiple outputs
    if force or not os.path.exists(os.path.join(core_dir, "unweighted_unifrac_pcoa_results.qza")):
        _run_qiime([qiime, "diversity", "core-metrics-phylogenetic",
              "--i-phylogeny", rooted_tree_qza,
              "--i-table", table_qza,
              "--p-sampling-depth", str(sampling_depth),
              "--output-dir", core_dir], log_file=log_path)

    _provenance(out, {"stage": "phylogeny_and_core_diversity", "sampling_depth": sampling_depth})

    return {
        "aligned_rep_seqs_qza": os.path.abspath(aligned_qza),
        "masked_alignment_qza": os.path.abspath(masked_qza),
        "unrooted_tree_qza": os.path.abspath(tree_qza),
        "rooted_tree_qza": os.path.abspath(rooted_tree_qza),
        "diversity_dir": os.path.abspath(core_dir),
    }

# 5) 统计与导出
def q2_stats_and_export(
    table_qza: str,
    taxonomy_qza: str,
    output_dir: str,
    qiime_bin: str,
    metadata_tsv: Optional[str] = None,
    methods: Tuple[str, ...] = ("ancom",),  # e.g., ('ancom', 'aldex2', 'songbird')
    export_biom: bool = True,
    export_tsv: bool = True,
    unpack_qzv: bool = True,
    force: bool = False,
) -> Dict[str, str]:
    """
    Run common differential abundance methods (if available), and export feature table and taxonomy
    to BIOM/TSV. Optionally unpack .qzv to static HTML assets for sharing.
    """
    table_qza = _require_file(table_qza, "table_qza")
    taxonomy_qza = _require_file(taxonomy_qza, "taxonomy_qza")
    qiime = _require_file(qiime_bin, "qiime_bin")
    if metadata_tsv:
        _require_file(metadata_tsv, "metadata_tsv")
    out = _ensure_dir(output_dir)
    log_path = os.path.join(out, "run.log")

    results: Dict[str, str] = {}

    # ANCOM (composition)
    if "ancom" in methods:
        comp_dir = os.path.join(out, "composition_ancom")
        _ensure_dir(comp_dir)
        comp_table_qza = os.path.join(comp_dir, "comp_table.qza")
        if force or not os.path.exists(comp_table_qza):
            _run_qiime([qiime, "composition", "add-pseudocount",
                  "--i-table", table_qza,
                  "--o-composition-table", comp_table_qza], log_file=log_path)
        if metadata_tsv:
            ancom_qzv = os.path.join(comp_dir, "ancom.qzv")
            # Require a metadata column name; we use a generic example 'Group' here,
            # caller should provide a metadata with a 'Group' column or adapt as needed.
            # To keep this wrapper generic, we record a note instead of guessing column.
            _provenance(comp_dir, {"note": "Run ANCOM via QIIME 2 manually specifying --m-metadata-column in downstream usage if needed."})
            results["ancom_composition_table_qza"] = os.path.abspath(comp_table_qza)
            results["ancom_qzv_placeholder"] = "Provide metadata column to run 'qiime composition ancom' manually."
        else:
            results["ancom_skipped"] = "No metadata_tsv provided."

    # ALDEx2 (if plugin installed)
    if "aldex2" in methods:
        ald_dir = os.path.join(out, "aldex2")
        _ensure_dir(ald_dir)
        _provenance(ald_dir, {"note": "Requires q2-aldex2 plugin; run specific commands outside or extend wrapper with column params."})
        results["aldex2_note"] = "q2-aldex2 not executed: column/conditions are dataset-specific."

    # Songbird (if plugin installed)
    if "songbird" in methods:
        sb_dir = os.path.join(out, "songbird")
        _ensure_dir(sb_dir)
        _provenance(sb_dir, {"note": "Requires q2-songbird; design formula and parameters are dataset-specific."})
        results["songbird_note"] = "q2-songbird not executed: design formula required."

    # Export BIOM/TSV
    if export_biom:
        export_biom_dir = os.path.join(out, "export_biom")
        _ensure_dir(export_biom_dir)
        if force or not os.path.exists(os.path.join(export_biom_dir, "feature-table.biom")):
            _run_qiime([qiime, "tools", "export",
                  "--input-path", table_qza,
                  "--output-path", export_biom_dir], log_file=log_path)
        biom_path = os.path.join(export_biom_dir, "feature-table.biom")
        results["biom_path"] = os.path.abspath(biom_path)

    if export_tsv:
        # Convert taxonomy to TSV
        tax_tsv_dir = os.path.join(out, "export_taxonomy")
        _ensure_dir(tax_tsv_dir)
        if force or not os.path.exists(os.path.join(tax_tsv_dir, "taxonomy.tsv")):
            _run_qiime([qiime, "tools", "export",
                  "--input-path", taxonomy_qza,
                  "--output-path", tax_tsv_dir], log_file=log_path)
        results["taxonomy_tsv"] = os.path.abspath(os.path.join(tax_tsv_dir, "taxonomy.tsv"))

        # Summarize feature table to a TSV using 'feature-table summarize' (qzv) + manual unpack is needed.
        # Here we export raw BIOM above; TSV conversion can be done with biom-format CLI (external),
        # so we record a note to keep stdlib-only constraint.
        _provenance(out, {"note": "For TSV feature table, convert BIOM via 'biom convert' (external) or Pandas in Python."})

    if unpack_qzv:
        _provenance(out, {"note": "To unpack .qzv visualizations, use 'qiime tools export' on each .qzv as needed."})

    _provenance(out, {"stage": "stats_and_export", "methods": methods})
    return results


def run_mash_deduplicate(
    fasta_paths: List[str],
    out_dir: str,
    mash_path: str = "./operation_env/tool_lake/mash",
    kmer: int = 21,
    sketch_size: int = 1000,
    max_distance: float = 0.004,
    threads: int = 1,
    keep_temps: bool = True,
    distance_tsv_name: str = "mash_self_dist.tsv",
) -> Dict[str, object]:
    """
    Run Mash-based de-replication on a set of genome FASTA files.

    This function builds a Mash sketch for the provided genomes and computes
    all-vs-all Mash distances, then clusters genomes by single-linkage under
    a distance threshold (default 0.004 ≈ 99.6% identity). One representative
    per cluster is selected (the largest file size within the cluster).

    Args:
        fasta_paths (List[str]): Paths to input genome FASTA files.
        mash_path (str): Absolute path to the Mash executable (e.g., "/opt/mash/mash").
                         The function does NOT rely on system PATH.
        out_dir (str): Output directory to store sketch, distance table, and logs.
        kmer (int, optional): Mash k-mer size for sketching. Defaults to 21.
        sketch_size (int, optional): Mash sketch size (-s). Defaults to 1000.
        max_distance (float, optional): Max Mash distance to consider two genomes
                                        redundant (clustered together). Defaults to 0.004
                                        (≈ 99.6% identity).
        threads (int, optional): Number of CPU threads for Mash where applicable.
                                 Defaults to 1.
        keep_temps (bool, optional): Keep intermediate files (.msh, temp list). Defaults to True.
        distance_tsv_name (str, optional): Filename for the Mash distance table (.tsv). Defaults to "mash_self_dist.tsv".

    Returns:
        Dict[str, object]: A dictionary with keys:
            - "representatives": List[str] of chosen representative FASTA paths.
            - "clusters": List[List[str]] where each sublist is the filenames in a cluster.
            - "distance_table": str, path to the written Mash distance table (TSV).
            - "sketch_path": str, path to the Mash sketch (.msh).
            - "params": Dict of parameters used (kmer, sketch_size, max_distance, threads).

    Raises:
        FileNotFoundError: If mash_path or any FASTA is not found.
        ValueError: If input list is empty or parameters are invalid.
        RuntimeError: If Mash invocation fails.

    Notes:
        - Mash 'distance' column is used. For highly similar genomes, 1 - distance
          approximates ANI. A 0.004 threshold roughly matches the 99.6% similarity cutoff
          used in一些数据库去冗余实践。
        - The procedure is O(N^2) in number of genomes because it computes all-vs-all distances.
          For very large N, consider pre-binning or sketch down-sampling.
    Examples:
        >>> result = run_mash_deduplicate(
        ...     fasta_paths=["A.fna", "B.fna", "C.fna"],
        ...     mash_path="/opt/mash/mash",
        ...     out_dir="./mash_out",
        ...     max_distance=0.004
        ... )
        >>> result["representatives"]
        ['B.fna', 'C.fna']
    """
    # -------------------- preflight checks --------------------
    if not isinstance(fasta_paths, list) or len(fasta_paths) == 0:
        raise ValueError("fasta_paths must be a non-empty list of FASTA file paths.")
    if not isinstance(mash_path, str) or not os.path.isfile(mash_path):
        raise FileNotFoundError(f"mash executable not found at: {mash_path}")
    for fp in fasta_paths:
        if not os.path.isfile(fp):
            raise FileNotFoundError(f"FASTA not found: {fp}")
    if kmer <= 0 or sketch_size <= 0 or max_distance < 0.0:
        raise ValueError("Invalid parameters: kmer/sketch_size must be >0, max_distance >= 0.")
    os.makedirs(out_dir, exist_ok=True)

    # -------------------- prepare file list --------------------
    # Mash can accept a list of files; we'll pass them directly to 'mash sketch'
    tmpdir = tempfile.mkdtemp(prefix="mash_dedupe_", dir=out_dir)
    list_path = os.path.join(tmpdir, "genomes.list")
    with open(list_path, "w", encoding="utf-8") as f:
        for fp in fasta_paths:
            f.write(os.path.abspath(fp) + "\n")

    sketch_path = os.path.join(out_dir, "genomes.msh")
    distance_tsv = os.path.join(out_dir, distance_tsv_name)
    log_path = os.path.join(out_dir, "mash.log")

    # -------------------- run mash sketch --------------------
    # mash sketch -k <kmer> -s <sketch_size> -o <sketch_path> -l <list_path>
    sketch_cmd = [
        mash_path, "sketch",
        "-k", str(kmer),
        "-s", str(sketch_size),
        "-o", sketch_path,
        "-l", list_path
    ]
    # Some mash builds accept -p for threads in distance, not always in sketch; safe to omit here.

    try:
        with open(log_path, "a", encoding="utf-8") as logf:
            logf.write("[CMD] " + " ".join(sketch_cmd) + "\n")
            proc = subprocess.run(sketch_cmd, stdout=logf, stderr=logf, check=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"mash sketch failed (see log: {log_path})") from e

    # Mash outputs sketch_path.msh; ensure extension
    if not os.path.isfile(sketch_path + ".msh"):
        # Some versions require no extra extension in -o; standard is to append ".msh"
        # Normalize the path to the actual .msh
        if os.path.isfile(sketch_path):
            msh_file = sketch_path
        else:
            # try to find a .msh in out_dir
            candidates = [p for p in os.listdir(out_dir) if p.endswith(".msh")]
            if not candidates:
                raise RuntimeError("Mash sketch did not produce a .msh file.")
            msh_file = os.path.join(out_dir, candidates[0])
    else:
        msh_file = sketch_path + ".msh"

    # -------------------- run mash dist (self vs self) --------------------
    # mash dist -p <threads> genomes.msh genomes.msh > distances.tsv
    dist_cmd = [mash_path, "dist", "-p", str(threads), msh_file, msh_file]
    try:
        with open(distance_tsv, "w", encoding="utf-8") as outf, \
             open(log_path, "a", encoding="utf-8") as logf:
            logf.write("[CMD] " + " ".join(dist_cmd) + f" > {distance_tsv}\n")
            proc = subprocess.run(dist_cmd, stdout=outf, stderr=logf, check=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"mash dist failed (see log: {log_path})") from e

    if not os.path.isfile(distance_tsv):
        raise RuntimeError("Expected distance TSV not found after mash dist.")

    # -------------------- parse distance table --------------------
    # mash dist output columns: query, ref, distance, p-value, shared-hashes
    # We'll build an adjacency list for pairs with distance <= max_distance (excluding self-pairs).
    index_by_path = {os.path.abspath(fp): i for i, fp in enumerate(fasta_paths)}
    n = len(fasta_paths)
    adj = [[] for _ in range(n)]

    def _norm(path: str) -> str:
        # Mash outputs the original path used in sketch; normalize to abs
        return os.path.abspath(path)

    with open(distance_tsv, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 3:
                # Some builds use spaces; try split by any whitespace
                parts = line.strip().split()
                if len(parts) < 3:
                    continue
            q, r, dist = parts[0], parts[1], parts[2]
            try:
                d = float(dist)
            except ValueError:
                continue
            qn, rn = _norm(q), _norm(r)
            if qn == rn:
                continue
            if qn in index_by_path and rn in index_by_path and d <= max_distance:
                qi, ri = index_by_path[qn], index_by_path[rn]
                adj[qi].append(ri)
                adj[ri].append(qi)

    # -------------------- connected components (single-linkage) --------------------
    visited = [False] * n
    clusters: List[List[str]] = []

    for i in range(n):
        if visited[i]:
            continue
        # BFS/DFS
        stack = [i]
        visited[i] = True
        comp_idx = [i]
        while stack:
            u = stack.pop()
            for v in adj[u]:
                if not visited[v]:
                    visited[v] = True
                    stack.append(v)
                    comp_idx.append(v)
        # convert to paths
        clusters.append([os.path.abspath(fasta_paths[j]) for j in comp_idx])

    # -------------------- choose representatives (largest file size per cluster) --------------------
    def _file_size(p: str) -> int:
        try:
            return os.path.getsize(p)
        except OSError:
            return -1

    representatives: List[str] = []
    for comp in clusters:
        rep = max(comp, key=_file_size)
        representatives.append(rep)

    # -------------------- cleanup --------------------
    if not keep_temps:
        try:
            shutil.rmtree(tmpdir, ignore_errors=True)
        except Exception:
            pass

    return {
        "representatives": representatives,
        "clusters": clusters,
        "distance_table": distance_tsv,
        "sketch_path": msh_file,
        "params": {
            "kmer": kmer,
            "sketch_size": sketch_size,
            "max_distance": max_distance,
            "threads": threads,
        },
    }


# 1) ARTS（抗生素耐药靶标搜寻器）
def run_arts(input_fasta: str, output_dir: str, arts_exec: str = "./operation_env/tool_lake/arts", extra_args: Optional[List[str]] = None) -> str:
    """
    调用 ARTS 命令行进行分析。
    Args:
        arts_exec: ARTS 可执行文件的路径（如 /opt/arts/bin/arts）
        input_fasta: 输入序列（fasta）
        output_dir: 输出目录
        extra_args: 追加的原生命令行参数列表
    Returns:
        输出目录路径（存在且包含结果）
    """
    _ensure_file(arts_exec, "ARTS executable")
    _ensure_file(input_fasta, "Input FASTA")
    _ensure_dir(output_dir, "Output directory")
    cmd = [arts_exec, "-i", input_fasta, "-o", output_dir]
    if extra_args:
        cmd += extra_args
    _run(cmd)
    return output_dir

# 2) Exonerate
def run_exonerate(query_fasta: str, target_fasta: str,
                  model: str, out_file: str, exonerate_path: str = "./operation_env/tool_lake/exonerate", extra_args: Optional[List[str]] = None) -> str:
    """
    运行 Exonerate 比对。
    """
    _ensure_file(exonerate_path, "Exonerate executable")
    _ensure_file(query_fasta, "Query FASTA")
    _ensure_file(target_fasta, "Target FASTA")
    out_dir = os.path.dirname(out_file) or "."
    _ensure_dir(out_dir, "Output dir")
    cmd = [exonerate_path, f"--model={model}", query_fasta, target_fasta]
    if extra_args:
        cmd += extra_args
    res = _run(cmd)
    with open(out_file, "w", encoding="utf-8") as f:
        f.write(res.stdout)
    return out_file

# 3) Cactus（多基因组比对）
def run_cactus(job_store: str, seqfile: str, out_hal: str, cactus_path: str = "./operation_env/tool_lake/cactus", 
               extra_args: Optional[List[str]] = None) -> str:
    """
    运行 Cactus 主流程。
    """
    _ensure_file(cactus_path, "Cactus executable")
    _ensure_file(seqfile, "Cactus seqfile")
    _ensure_dir(os.path.dirname(out_hal) or ".", "Output dir")
    cmd = [cactus_path, job_store, seqfile, out_hal]
    if extra_args:
        cmd += extra_args
    _run(cmd)
    if not os.path.isfile(out_hal):
        raise RuntimeError(f"Cactus expected output not found: {out_hal}")
    return out_hal

# 4) FastQC（测序质控）
def run_fastqc(inputs: List[str], out_dir: str, threads: int = 1, fastqc_path: str = "./operation_env/tool_lake/fastqc", 
               extra_args: Optional[List[str]] = None) -> str:
    """
    批量运行 FastQC。
    """
    _ensure_file(fastqc_path, "FastQC executable")
    for p in inputs:
        _ensure_file(p, f"Input file {p}")
    _ensure_dir(out_dir, "Output directory")
    cmd = [fastqc_path, "--threads", str(threads), "--outdir", out_dir] + inputs
    if extra_args:
        cmd += extra_args
    _run(cmd)
    return out_dir

# 5) Trimmomatic（Java JAR）
def run_trimmomatic(mode: str,
                    input_R1: str, input_R2: Optional[str],
                    output_R1: str, output_R1_unpaired: Optional[str],
                    output_R2: Optional[str], output_R2_unpaired: Optional[str],
                    trimmomatic_jar: str = "./operation_env/tool_lake/trimmomatic",
                    adapters_fa: Optional[str] = None,
                    threads: int = 4, extra_args: Optional[List[str]] = None, java_path: str = "java",) -> List[str]:
    """
    运行 Trimmomatic（支持 PE/SE）。
    mode: "PE" 或 "SE"
    """
    # _ensure_file(java_path, "Java executable")
    _ensure_file(trimmomatic_jar, "Trimmomatic JAR")
    _ensure_file(input_R1, "Input R1")
    if mode.upper() == "PE":
        if input_R2 is None:
            raise ValueError("PE mode requires input_R2")
        _ensure_file(input_R2, "Input R2")
        outputs = [output_R1, output_R1_unpaired, output_R2, output_R2_unpaired]
        if any(o is None for o in outputs):
            raise ValueError("PE mode requires all four output paths")
        cmd = [java_path, "-jar", trimmomatic_jar, "PE", "-threads", str(threads),
               input_R1, input_R2, output_R1, output_R1_unpaired, output_R2, output_R2_unpaired]
    elif mode.upper() == "SE":
        outputs = [output_R1]
        cmd = [java_path, "-jar", trimmomatic_jar, "SE", "-threads", str(threads),
               input_R1, output_R1]
    else:
        raise ValueError("mode must be 'PE' or 'SE'")
    if adapters_fa:
        _ensure_file(adapters_fa, "Adapters FASTA")
        cmd += ["ILLUMINACLIP:" + adapters_fa + ":2:30:10"]
    if extra_args:
        cmd += extra_args
    _run(cmd)
    return [o for o in outputs if o]

# 6) NCBI datasets（命令行）
def run_ncbi_datasets(subcommand: List[str], out_dir: str, datasets_exec: str = "./operation_env/tool_lake/ncbidatasets", ) -> str:
    """
    通用封装：例如下载基因组
    subcommand 例：["download", "genome", "accession", "--input-file", "acc.txt", "--filename", "dl.zip"]
    """
    _ensure_file(datasets_exec, "datasets executable")
    _ensure_dir(out_dir, "Output directory")
    cmd = [datasets_exec] + subcommand
    _run(cmd, cwd=out_dir)
    return out_dir

# 7) NCBI dataformat（命令行）
def run_ncbi_dataformat(subcommand: List[str], dataformat_exec: str = "./operation_env/tool_lake/dataformat", workdir: Optional[str] = None) -> subprocess.CompletedProcess:
    """
    通用封装：将 datasets 的 JSONL 转 TSV 等
    subcommand 例：["tsv", "genome", "--fields", "accession,organism-name", "--inputfile", "data.jsonl"]
    """
    _ensure_file(dataformat_exec, "dataformat executable")
    return _run([dataformat_exec] + subcommand, cwd=workdir)

# 8) InterProScan
def run_interproscan(input_fasta: str, out_dir: str,
                     formats: str = "TSV", interproscan_exec: str = "./operation_env/tool_lake/InterProScan",  applications: Optional[List[str]] = None,
                     extra_args: Optional[List[str]] = None) -> str:
    """
    运行 InterProScan。
    """
    _ensure_file(interproscan_exec, "InterProScan executable")
    _ensure_file(input_fasta, "Input FASTA")
    _ensure_dir(out_dir, "Output directory")
    cmd = [interproscan_exec, "-i", input_fasta, "-f", formats, "-d", out_dir]
    if applications:
        cmd += ["-appl", ",".join(applications)]
    if extra_args:
        cmd += extra_args
    _run(cmd)
    return out_dir

# 10) Bowtie2（比对）
def run_bowtie2(index_base: str, read1: str,
                out_sam: str, bowtie2_exec: str = "./operation_env/tool_lake/Bowtie2", read2: Optional[str] = None, threads: int = 4,
                extra_args: Optional[List[str]] = None) -> str:
    """
    运行 Bowtie2，支持 SE/PE。
    """
    _ensure_file(bowtie2_exec, "bowtie2 executable")
    _ensure_file(read1, "Read1")
    if read2:
        _ensure_file(read2, "Read2")
    _ensure_dir(os.path.dirname(out_sam) or ".", "Output dir")
    cmd = [bowtie2_exec, "-x", index_base, "-p", str(threads)]
    if read2:
        cmd += ["-1", read1, "-2", read2]
    else:
        cmd += ["-U", read1]
    if extra_args:
        cmd += extra_args
    res = _run(cmd)
    with open(out_sam, "w", encoding="utf-8") as f:
        f.write(res.stdout)
    return out_sam


# 12) ANNOVAR（table_annovar.pl）
def run_annovar(table_annovar_pl: str, annovar_dir: str, humandb_dir: str,
                input_vcf: str, buildver: str, out_prefix: str,
                protocols: List[str], operations: List[str],
                other_args: Optional[List[str]] = None) -> List[str]:
    """
    运行 ANNOVAR 注释（table_annovar.pl）。
    """
    _ensure_file(table_annovar_pl, "table_annovar.pl")
    _ensure_dir(annovar_dir, "annovar_dir")
    _ensure_dir(humandb_dir, "humandb_dir")
    _ensure_file(input_vcf, "Input VCF")
    out_dir = os.path.dirname(out_prefix) or "."
    _ensure_dir(out_dir, "Output dir")

    cmd = [table_annovar_pl, input_vcf, humandb_dir,
           "-buildver", buildver,
           "-out", out_prefix,
           "-remove",
           "-protocol", ",".join(protocols),
           "-operation", ",".join(operations),
           "-nastring", ".",
           "-vcfinput"]
    # 让脚本能找到 annovar 自带 perl 脚本/数据库（有些环境需要）
    env = os.environ.copy()
    env["ANNOVAR_DIR"] = annovar_dir
    if other_args:
        cmd += other_args

    res = subprocess.run(cmd, capture_output=True, text=True, check=False, env=env)
    if res.returncode != 0:
        raise RuntimeError(
            "ANNOVAR failed:\n"
            f"$ {' '.join(shlex.quote(c) for c in cmd)}\n"
            f"STDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}"
        )

    # 简单枚举常见输出：multianno.* 或 *_hg19_multianno.txt 等
    generated = []
    for fn in os.listdir(out_dir):
        if os.path.basename(out_prefix) in fn:
            generated.append(os.path.join(out_dir, fn))
    return sorted(generated)