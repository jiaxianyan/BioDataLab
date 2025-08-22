import os
import json
import shlex
import subprocess
import datetime
from typing import Optional, Dict, Any, List, Tuple, Sequence, Union
import tempfile
import shutil

import qiime2

# QIIME 2 插件
from qiime2.plugins import (
    tools,
    demux,
    vsearch,
    cutadapt,
    quality_filter,
    deblur,
    dada2,
    feature_table,
    feature_classifier,
    taxa,
    alignment,
    phylogeny,
    diversity,
    c
)

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
def q2_ingest_and_qc_api(
    manifest_csv: str,
    output_dir: str,
    sample_metadata_tsv: str = None,
    paired_end: bool = True,
    join_paired: bool = False,
    cutadapt_params: dict = None,   # 改成 dict 传参，例如 {"p_front_f": "ACTG"}
    qscore_filter: bool = False,
    trunc_len_f: int = None,
    trunc_len_r: int = None,
    phred_offset: int = 33,
    force: bool = False,
):
    """
    Python API 版：导入 FASTQ（manifest），可选 join-pairs, cutadapt 修剪, q-score 过滤。
    返回 demux.qza 和 demux.qzv 的 Artifact 对象。
    """
    os.makedirs(output_dir, exist_ok=True)

    fmt = ("PairedEndFastqManifestPhred{}V2" if paired_end else "SingleEndFastqManifestPhred{}V2").format(phred_offset)

    # 导入数据
    demux_art = tools.import_(
        type="SampleData[PairedEndSequencesWithQuality]" if paired_end else "SampleData[SequencesWithQuality]",
        input_path=manifest_csv,
        input_format=fmt
    ).result

    work_art = demux_art

    # join-pairs
    if paired_end and join_paired:
        joined = vsearch.join_pairs(i_demultiplexed_seqs=work_art)
        work_art = joined.o_joined_sequences

    # cutadapt 修剪
    if cutadapt_params:
        if paired_end:
            trimmed = cutadapt.trim_paired(i_demultiplexed_sequences=work_art, **cutadapt_params)
        else:
            trimmed = cutadapt.trim_single(i_demultiplexed_sequences=work_art, **cutadapt_params)
        work_art = trimmed.o_trimmed_sequences

    # q-score filter
    if qscore_filter:
        if paired_end and join_paired:
            filt = quality_filter.q_score_joined(i_demux=work_art)
        else:
            filt = quality_filter.q_score(i_demux=work_art)
        work_art = filt.o_filtered_sequences

    # summarize
    summary = demux.summarize(i_data=work_art)

    # 保存输出
    demux_qza = os.path.join(output_dir, "demux.qza")
    demux_qzv = os.path.join(output_dir, "demux.qzv")
    work_art.save(demux_qza)
    summary.visualization.save(demux_qzv)

    return {"demux_qza": demux_qza, "demux_qzv": demux_qzv}

# 2) 去噪与特征表
def q2_denoise_and_feature_table_api(
    demux_qza: str,
    output_dir: str,
    method: str = "deblur",   # 'deblur' or 'dada2'
    trim_length: int = 250,   # for deblur
    trunc_len_f: int = None,  # for dada2 paired
    trunc_len_r: int = None,  # for dada2 paired
    min_reads_per_sample: int = None,
    threads: int = 0,
    force: bool = False,
):
    """
    Python API 版：用 Deblur 或 DADA2 去噪，生成特征表和代表序列，
    可选过滤低丰度样本。
    """
    os.makedirs(output_dir, exist_ok=True)

    # 加载输入 artifact
    demux_art = qiime2.Artifact.load(demux_qza)

    # 输出文件
    table_qza = os.path.join(output_dir, "table.qza")
    rep_qza   = os.path.join(output_dir, "rep-seqs.qza")
    stats_qzv = os.path.join(output_dir, f"{method}_stats.qzv")

    # ------------------- 去噪 -------------------
    if method.lower() == "deblur":
        result = deblur.denoise_16S(
            i_demultiplexed_seqs=demux_art,
            p_trim_length=trim_length,
            p_jobs_to_start=threads if threads > 0 else None
        )
        rep_seqs = result.o_representative_sequences
        table    = result.o_table
        stats    = result.o_stats

    elif method.lower() == "dada2":
        if trunc_len_f is None and trunc_len_r is None:
            # single-end
            result = dada2.denoise_single(
                i_demultiplexed_seqs=demux_art,
                p_trunc_len=trunc_len_f if trunc_len_f else 0,
                p_n_threads=threads if threads > 0 else 0
            )
        else:
            # paired-end
            result = dada2.denoise_paired(
                i_demultiplexed_seqs=demux_art,
                p_trunc_len_f=trunc_len_f if trunc_len_f else 0,
                p_trunc_len_r=trunc_len_r if trunc_len_r else 0,
                p_n_threads=threads if threads > 0 else 0
            )
        rep_seqs = result.o_representative_sequences
        table    = result.o_table
        stats    = result.o_denoising_stats

    else:
        raise ValueError("method must be 'deblur' or 'dada2'.")

    # 保存主要结果
    rep_seqs.save(rep_qza)
    table.save(table_qza)
    stats.visualization.save(stats_qzv)

    # ------------------- 可选：过滤低丰度样本 -------------------
    filtered_table_qza = ""
    if min_reads_per_sample and min_reads_per_sample > 0:
        filtered_table = feature_table.filter_samples(
            i_table=table,
            p_min_frequency=min_reads_per_sample
        ).o_filtered_table
        filtered_table_qza = os.path.join(output_dir, f"table_min{min_reads_per_sample}.qza")
        filtered_table.save(filtered_table_qza)

    # ------------------- 返回结果 -------------------
    result_paths = {
        "table_qza": os.path.abspath(filtered_table_qza or table_qza),
        "rep_seqs_qza": os.path.abspath(rep_qza),
        "denoise_stats_qzv": os.path.abspath(stats_qzv),
    }
    if filtered_table_qza:
        result_paths["table_filtered_qza"] = os.path.abspath(filtered_table_qza)

    return result_paths

# 3) 分类学注释
def q2_taxonomy_api(
    table_qza: str,
    rep_seqs_qza: str,
    output_dir: str,
    metadata_tsv: str = None,
    pretrained_classifier_qza: str = None,
    ref_seqs_qza: str = None,
    ref_taxonomy_qza: str = None,
    method: str = "sklearn",   # 'sklearn' or 'vsearch'
    perc_identity: float = 0.97,
    force: bool = False,
):
    """
    Python API 版：使用预训练分类器或参考数据库进行分类学注释，
    然后生成 taxa barplot。
    """
    os.makedirs(output_dir, exist_ok=True)

    # 加载输入
    table_art = qiime2.Artifact.load(table_qza)
    rep_seqs_art = qiime2.Artifact.load(rep_seqs_qza)
    metadata = qiime2.Metadata.load(metadata_tsv) if metadata_tsv else None

    taxonomy_qza = os.path.join(output_dir, "taxonomy.qza")
    taxonomy_qzv = os.path.join(output_dir, "taxonomy.qzv")
    taxa_barplot_qzv = os.path.join(output_dir, "taxa-barplot.qzv")

    classifier_art = None

    # ------------------- 选择分类器 -------------------
    if pretrained_classifier_qza:
        classifier_art = qiime2.Artifact.load(pretrained_classifier_qza)
    else:
        if not (ref_seqs_qza and ref_taxonomy_qza):
            raise ValueError("Either provide 'pretrained_classifier_qza' OR both 'ref_seqs_qza' and 'ref_taxonomy_qza'.")
        ref_seqs_art = qiime2.Artifact.load(ref_seqs_qza)
        ref_tax_art  = qiime2.Artifact.load(ref_taxonomy_qza)

        if method == "sklearn":
            # 训练 Naive Bayes 分类器
            classifier_art = feature_classifier.fit_classifier_naive_bayes(
                i_reference_reads=ref_seqs_art,
                i_reference_taxonomy=ref_tax_art
            ).o_classifier
            classifier_art.save(os.path.join(output_dir, "trained_classifier.qza"))

    # ------------------- 分类学注释 -------------------
    if method == "sklearn":
        result = feature_classifier.classify_sklearn(
            i_classifier=classifier_art,
            i_reads=rep_seqs_art
        )
        taxonomy = result.o_classification

    elif method == "vsearch":
        if not (ref_seqs_qza and ref_taxonomy_qza):
            raise ValueError("VSEARCH method requires 'ref_seqs_qza' and 'ref_taxonomy_qza'.")
        ref_seqs_art = qiime2.Artifact.load(ref_seqs_qza)
        ref_tax_art  = qiime2.Artifact.load(ref_taxonomy_qza)

        result = feature_classifier.classify_consensus_vsearch(
            i_query=rep_seqs_art,
            i_reference_reads=ref_seqs_art,
            i_reference_taxonomy=ref_tax_art,
            p_perc_identity=perc_identity
        )
        taxonomy = result.o_classification

    else:
        raise ValueError("method must be 'sklearn' or 'vsearch'.")

    taxonomy.save(taxonomy_qza)

    # ------------------- 生成 barplot -------------------
    barplot = taxa.barplot(
        i_table=table_art,
        i_taxonomy=taxonomy,
        m_metadata_file=metadata
    )
    barplot.visualization.save(taxa_barplot_qzv)

    # 保存 taxonomy 可视化
    taxonomy.visualization = feature_classifier.classify_sklearn(  # trick: reuse for viz
        i_classifier=classifier_art, i_reads=rep_seqs_art
    ).o_classification.view(qiime2.Metadata)

    # ------------------- 返回结果 -------------------
    result_paths = {
        "taxonomy_qza": os.path.abspath(taxonomy_qza),
        "taxonomy_qzv": os.path.abspath(taxonomy_qzv),
        "taxa_barplot_qzv": os.path.abspath(taxa_barplot_qzv),
    }
    if pretrained_classifier_qza is None and method == "sklearn":
        result_paths["trained_classifier_qza"] = os.path.abspath(
            os.path.join(output_dir, "trained_classifier.qza")
        )

    return result_paths

# 4) 系统发育与核心多样性
def q2_phylogeny_and_core_diversity_api(
    table_qza: str,
    rep_seqs_qza: str,
    output_dir: str,
    sampling_depth: int = 1000,
    use_midpoint_root: bool = True,
    force: bool = False,
):
    """
    Python API 版：构建系统发育树，并运行 core-metrics-phylogenetic。
    返回比对、树和多样性分析结果路径。
    """
    os.makedirs(output_dir, exist_ok=True)

    # 加载输入
    table_art = qiime2.Artifact.load(table_qza)
    rep_seqs_art = qiime2.Artifact.load(rep_seqs_qza)

    # 输出文件路径
    aligned_qza     = os.path.join(output_dir, "aligned-rep-seqs.qza")
    masked_qza      = os.path.join(output_dir, "masked-alignment.qza")
    tree_qza        = os.path.join(output_dir, "unrooted-tree.qza")
    rooted_tree_qza = os.path.join(output_dir, "rooted-tree.qza")
    core_dir        = os.path.join(output_dir, "core-metrics")
    os.makedirs(core_dir, exist_ok=True)

    # ------------------- 比对 -------------------
    aligned = alignment.mafft(i_sequences=rep_seqs_art).o_alignment
    aligned.save(aligned_qza)

    # ------------------- 掩码 -------------------
    masked = alignment.mask(i_alignment=aligned).o_masked_alignment
    masked.save(masked_qza)

    # ------------------- 构树 -------------------
    tree = phylogeny.fasttree(i_alignment=masked).o_tree
    tree.save(tree_qza)

    # ------------------- 加根 -------------------
    if use_midpoint_root:
        rooted = phylogeny.midpoint_root(i_tree=tree).o_rooted_tree
    else:
        rooted = phylogeny.root(i_tree=tree).o_rooted_tree
    rooted.save(rooted_tree_qza)

    # ------------------- 核心多样性指标 -------------------
    core = diversity.core_metrics_phylogenetic(
        i_table=table_art,
        i_phylogeny=rooted,
        p_sampling_depth=sampling_depth,
        output_dir=core_dir
    )

    # ------------------- 返回结果 -------------------
    return {
        "aligned_rep_seqs_qza": os.path.abspath(aligned_qza),
        "masked_alignment_qza": os.path.abspath(masked_qza),
        "unrooted_tree_qza": os.path.abspath(tree_qza),
        "rooted_tree_qza": os.path.abspath(rooted_tree_qza),
        "diversity_dir": os.path.abspath(core_dir),
    }

# 5) 统计与导出
def q2_stats_and_export_api(
    table_qza: str,
    taxonomy_qza: str,
    output_dir: str,
    metadata_tsv: str = None,
    methods: tuple = ("ancom",),  # e.g., ('ancom', 'aldex2', 'songbird')
    export_biom: bool = True,
    export_tsv: bool = True,
    unpack_qzv: bool = True,
    force: bool = False,
):
    """
    Python API 版：执行常见差异丰度分析（如 ANCOM），并导出 feature table 与 taxonomy。
    返回结果文件路径。
    """
    os.makedirs(output_dir, exist_ok=True)

    # 加载输入
    table_art = qiime2.Artifact.load(table_qza)
    taxonomy_art = qiime2.Artifact.load(taxonomy_qza)
    metadata = qiime2.Metadata.load(metadata_tsv) if metadata_tsv else None

    results = {}

    # ------------------- ANCOM -------------------
    if "ancom" in methods:
        comp_dir = os.path.join(output_dir, "composition_ancom")
        os.makedirs(comp_dir, exist_ok=True)

        comp_table_qza = os.path.join(comp_dir, "comp_table.qza")
        if force or not os.path.exists(comp_table_qza):
            comp_table = composition.add_pseudocount(i_table=table_art).o_composition_table
            comp_table.save(comp_table_qza)

        results["ancom_composition_table_qza"] = os.path.abspath(comp_table_qza)
        results["ancom_note"] = (
            "要运行 ANCOM 需指定 metadata column，请在调用 composition.ancom 时传入 --m-metadata-column"
        )

    # ------------------- ALDEx2 占位 -------------------
    if "aldex2" in methods:
        ald_dir = os.path.join(output_dir, "aldex2")
        os.makedirs(ald_dir, exist_ok=True)
        results["aldex2_note"] = "需要 q2-aldex2 插件；请根据数据集具体设计参数。"

    # ------------------- Songbird 占位 -------------------
    if "songbird" in methods:
        sb_dir = os.path.join(output_dir, "songbird")
        os.makedirs(sb_dir, exist_ok=True)
        results["songbird_note"] = "需要 q2-songbird 插件；请根据数据集具体设计公式。"

    # ------------------- 导出 BIOM -------------------
    if export_biom:
        export_biom_dir = os.path.join(output_dir, "export_biom")
        os.makedirs(export_biom_dir, exist_ok=True)
        tools.export_(input=table_art, output_path=export_biom_dir)
        results["biom_path"] = os.path.abspath(os.path.join(export_biom_dir, "feature-table.biom"))

    # ------------------- 导出 taxonomy TSV -------------------
    if export_tsv:
        export_tax_dir = os.path.join(output_dir, "export_taxonomy")
        os.makedirs(export_tax_dir, exist_ok=True)
        tools.export_(input=taxonomy_art, output_path=export_tax_dir)
        results["taxonomy_tsv"] = os.path.abspath(os.path.join(export_tax_dir, "taxonomy.tsv"))

    # ------------------- unpack .qzv 占位 -------------------
    if unpack_qzv:
        results["unpack_note"] = "要解包 .qzv 可用 tools.export_ 输入 visualization artifact。"

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
                  model: str, out_file: str,
                  exonerate_exe: str = "exonerate",
                  extra_args: Optional[List[str]] = None) -> str:
    """
    运行 Exonerate 比对。

    Args:
        query_fasta (str): 查询序列 FASTA 文件路径。
        target_fasta (str): 目标序列 FASTA 文件路径。
        model (str): 比对模型 (例如 "protein2genome", "dna2dna")。
        out_file (str): 输出结果文件路径。
        exonerate_exe (str): Exonerate 可执行文件 (默认 'exonerate'，即 conda 安装环境中)。
        extra_args (List[str], optional): 额外命令行参数，例如 ["--showalignment", "yes"]。

    Returns:
        str: 结果文件路径。

    Raises:
        FileNotFoundError: 如果 query 或 target 文件不存在。
        RuntimeError: 如果 Exonerate 执行失败。
    """
    if not os.path.isfile(query_fasta):
        raise FileNotFoundError(f"Query FASTA not found: {query_fasta}")
    if not os.path.isfile(target_fasta):
        raise FileNotFoundError(f"Target FASTA not found: {target_fasta}")

    os.makedirs(os.path.dirname(out_file) or ".", exist_ok=True)

    cmd = [exonerate_exe, f"--model={model}", query_fasta, target_fasta]
    if extra_args:
        cmd.extend(extra_args)

    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"Exonerate failed with exit code {proc.returncode}\n"
            f"CMD: {' '.join(shlex.quote(c) for c in cmd)}\n"
            f"STDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
        )

    with open(out_file, "w", encoding="utf-8") as f:
        f.write(proc.stdout)

    return os.path.abspath(out_file)

# 3) Cactus（多基因组比对）
def run_cactus(job_store: str, seqfile: str, out_hal: str,
               cactus_exe: str = "cactus",  # 默认直接用 conda 环境里的 cactus
               extra_args: Optional[List[str]] = None) -> str:
    """
    运行 Cactus 主流程。

    Args:
        job_store (str): Cactus 作业存储目录（必须是新目录或可重启目录）。
        seqfile (str): 输入的序列描述文件（seqFile.txt）。
        out_hal (str): 输出 HAL 文件路径。
        cactus_exe (str): Cactus 可执行文件名或路径（默认使用 conda 安装的 "cactus"）。
        extra_args (List[str], optional): 额外命令行参数，例如 ["--maxCores", "8"]。

    Returns:
        str: 生成的 HAL 文件的绝对路径。

    Raises:
        FileNotFoundError: 如果 seqfile 不存在。
        RuntimeError: 如果 Cactus 运行失败或没有生成 HAL 文件。
    """
    if not os.path.isfile(seqfile):
        raise FileNotFoundError(f"Seqfile not found: {seqfile}")
    os.makedirs(os.path.dirname(out_hal) or ".", exist_ok=True)

    cmd = [cactus_exe, job_store, seqfile, out_hal]
    if extra_args:
        cmd.extend(extra_args)

    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"Cactus failed with exit code {proc.returncode}\n"
            f"CMD: {' '.join(shlex.quote(c) for c in cmd)}\n"
            f"STDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
        )

    if not os.path.isfile(out_hal):
        raise RuntimeError(f"Cactus expected output not found: {out_hal}")

    return os.path.abspath(out_hal)

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
def run_ncbi_datasets(subcommand: List[str], out_dir: str,
                      datasets_exe: str = "datasets") -> str:
    """
    Run NCBI Datasets CLI command.

    Args:
        subcommand (List[str]): 子命令，例如:
            ["download", "genome", "accession", "--input-file", "acc.txt", "--filename", "dl.zip"]
        out_dir (str): 输出目录。
        datasets_exe (str): 可执行程序名，默认 "datasets"（来自 conda 安装）。

    Returns:
        str: 输出目录的绝对路径。

    Raises:
        RuntimeError: 如果执行失败。
    """
    os.makedirs(out_dir, exist_ok=True)

    cmd = [datasets_exe] + subcommand
    proc = subprocess.run(cmd, cwd=out_dir, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"NCBI datasets failed with exit code {proc.returncode}\n"
            f"CMD: {' '.join(shlex.quote(c) for c in cmd)}\n"
            f"STDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
        )

    return os.path.abspath(out_dir)


def run_ncbi_dataformat(subcommand: List[str],
                        dataformat_exe: str = "dataformat",
                        workdir: Optional[str] = None) -> subprocess.CompletedProcess:
    """
    Run NCBI Dataformat CLI command.

    Args:
        subcommand (List[str]): 子命令，例如:
            ["tsv", "genome", "--fields", "accession,organism-name", "--inputfile", "data.jsonl"]
        dataformat_exe (str): 可执行程序名，默认 "dataformat"（来自 conda 安装）。
        workdir (str, optional): 运行目录。

    Returns:
        subprocess.CompletedProcess: 进程结果对象，包含 stdout/stderr。
    """
    cmd = [dataformat_exe] + subcommand
    proc = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"NCBI dataformat failed with exit code {proc.returncode}\n"
            f"CMD: {' '.join(shlex.quote(c) for c in cmd)}\n"
            f"STDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
        )
    return proc

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