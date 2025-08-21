太好了！下面我把这 **15 个更加多样化的一步型 benchmark 任务**，分别给出：

* 执行步骤（面向 agent 的可操作流程）
* 需要构建的“工具”（agent可调用的最小原子工具；给出命令模板/参数规范）
* 参考答案（refdata）如何准备（含版本锁定、抽样与度量）

为便于落地，我还先给出**统一目录与环境规范**，保证复现实验和自动评测。

---

# 统一规范（强烈建议照此组织）

**目录结构（示例）**

```
biodatalab/
  data_lake/           # 原始数据/外部下载（只追加不改写）
  work/                # 中间结果（可重算）
  benchmark/
    tasks/             # 任务定义（yaml/json）
    results/           # 运行产物（按任务索引/日期分目录）
    refdata/           # 参考答案（冻结）
  tools/               # agent工具封装脚本（python/cli）
  eval/                # 评测脚本
  meta/                # 版本锁定信息（工具版本、数据库版本、快照时间）
```

**环境（建议）**

* 使用 Conda + Mamba 管理环境；每个外部重型工具单独环境或容器（如 `gtdbtk`, `interproscan`, `exonerate`）。
* 每个任务记录：

  * 工具版本（如 `interproscan-5.65-97.0`, `gtdbtk-2.4.0`）
  * 数据库/知识库版本（如 `GTDB R220`, `KEGG YYYY-MM-DD 快照`）
  * 数据来源快照日期（GEO/SRA/BioProject 搜索日期）
* 参考答案（refdata）在首次跑通后**冻结**；保存 metadata（行数、MD5、工具和库版本）。

**评测通用策略**

* 结构化输出（TSV/JSON/FASTA 等），**排序+去重**后比对。
* 下载/检索类：比对 **清单（manifest）** 或 ID 列表一致性。
* 过滤/标准化类：比对 **元素集合** 或 **字段映射** 一致性。
* 注释/比对/富集类：比对 **表格键列** + **关键字段**（如 domain ID、taxonomy、路径名、align坐标）。
* 指标：`exact match`、`F1`、`Jaccard`、`行数偏差阈值` 等（每任务下文给出推荐）。

---

# Agent 工具库（可复用的 15 个原子工具）

为保持一步原子化，每个工具建议有统一接口：

```
run(tool_name, **kwargs) -> output_path_or_object
```

并在 `tools/` 下实现对应脚本（Python 封装 CLI）。下面每个任务会指明需要哪个工具与关键参数。你可将这些封装为：

* `ncbi_datasets_download_genome`
* `geo_download_matrix`
* `pubmed_fetch_abstracts`
* `seq_length_filter`
* `fastqc_trimmomatic`
* `idmap_to_hgnc`
* `pdb_mmcif_to_fasta_json`
* `gtdbtk_classify`
* `interproscan_annotate`
* `kegg_enrich`
* `arts_detect_bgc`
* `exonerate_protein2genome`
* `drugcentral_map_targets`
* `host_trait_join`
* `bioproject_keyword_search`

每个工具输出**单一产物**（或单一清单），便于评测。

---

# 15 个任务：步骤 | 工具 | 参考答案

## A. 数据下载 / 抽取类

### 1) 下载蓝细菌基因组（NCBI Assembly + `datasets`）

* **步骤**

  1. 用 `ncbi-datasets`（命令或API）按分类单元 *Cyanobacteria* 下载基因组 FASTA。
  2. 导出一个 `manifest.tsv`（包含 accession、物种名、文件路径、MD5）。
* **工具**：`ncbi_datasets_download_genome(taxon="Cyanobacteria", format="fasta")`

  * 命令模板：`datasets download genome taxon cyanobacteria --include genome --filename out.zip`
* **Refdata**

  * 使用固定日期快照（写入 `meta/snapshot.json`），解压后生成 **标准化清单**（排序、去重）。
  * 评测：`manifest` 全匹配（行数、accession 列、MD5）。

---

### 2) 下载 COVID-19 表达矩阵（GEO）

* **步骤**

  1. 用 `GEOquery`（R）或 `pydataverse`/自写脚本获取指定 GEO 系列（给出 ID 清单或关键词 `SARS-CoV-2` +人工确定 N 个代表数据集）。
  2. 导出统一格式的表达矩阵（基因为行，样本为列，保存 TSV）。
* **工具**：`geo_download_matrix(geo_ids=[...])`

  * R 模板：`getGEO("GSEXXXXX"); exprs(eset) -> tsv`
* **Refdata**

  * 锁定 GEO IDs 列表；保存标准矩阵（排序标准化）。
  * 评测：矩阵维度一致 + 行/列名集合完全一致 + 随机抽查5行数值差异 < 1e-6。

---

### 3) PubMed 摘要抽取（关键词检索）

* **步骤**

  1. `Entrez.esearch` 用 `("SARS-CoV-2" AND "immune escape")` 检索 PMIDs（限定日期区间）。
  2. `Entrez.efetch` 拉取摘要，保存 `pmid\t标题\t摘要`。
* **工具**：`pubmed_fetch_abstracts(query, date_from, date_to)`
* **Refdata**

  * 冻结 PMIDs 列表；保存抽取表。
  * 评测：PMIDs 集合一致；若 PubMed 更新可允许±N条浮动，使用 `Jaccard≥0.95`。

---

## B. 数据过滤 / 质量控制类

### 4) 蛋白序列长度过滤（RefSeq 50–500 aa）

* **步骤**

  1. 输入蛋白 FASTA；保留长度 `[50,500]` 的序列；输出过滤后的 FASTA。
  2. 同时产出一个 `ids.txt`（序列ID列表）。
* **工具**：`seq_length_filter(input_fasta, min=50, max=500)`
* **Refdata**

  * 对指定输入 FASTA 过滤一次，冻结输出 FASTA + `ids.txt`。
  * 评测：`ids.txt` 全匹配（顺序无关），FASTA 记录数一致。

---

### 5) 宏基因组 reads 质控（FastQC + Trimmomatic）

* **步骤**

  1. `fastqc` 生成质控报告（可选评测）。
  2. `Trimmomatic` 按固定参数裁剪（如 `SLIDINGWINDOW:4:20 MINLEN:50`），输出 clean FASTQ。
* **工具**：`fastqc_trimmomatic(input_fastq, params=...)`
* **Refdata**

  * 小样本（如 50k reads）跑一次，冻结 clean FASTQ（或 `ids_kept.txt`）。
  * 评测：保留 reads 数一致；可比较 `read count` 和首尾 3 条 reads 的 MD5。

---

## C. 数据标准化 / 格式转换类

### 6) 基因 ID → HGNC 标准化（GEO 矩阵）

* **步骤**

  1. 输入表达矩阵（行名为多种基因 ID）。
  2. 用映射表（如 `org.Hs.eg.db` 或 Ensembl 映射）统一到 HGNC Symbol；丢弃无法映射条目；按 Symbol 聚合（mean/median）。
* **工具**：`idmap_to_hgnc(expr_tsv, map_db="org.Hs.eg.db")`
* **Refdata**

  * 冻结映射版本（包版本+日期）；冻结标准化矩阵。
  * 评测：行名集合完全一致；聚合策略固定（写入 meta）；随机抽查5个基因聚合值一致（±1e-6）。

---

### 7) mmCIF → FASTA + JSON（PDB 结构标准化）

* **步骤**

  1. 输入一组 PDB mmCIF 结构。
  2. 提取链级别氨基酸序列为 FASTA；将金属离子/配体等元信息整理为 JSON。
* **工具**：`pdb_mmcif_to_fasta_json(input_dir_or_list)`

  * Python：`Biopython` 的 `MMCIFParser`、`PPBuilder`。
* **Refdata**

  * 固定一组 PDB IDs；冻结 FASTA 与 JSON（键顺序与缩进固定）。
  * 评测：FASTA 条目集完全一致；JSON 用 `jsondiff` 比较（忽略字段顺序）。

---

### 8) GTDB-Tk 分类标准化（细菌基因组）

* **步骤**

  1. 输入一组基因组 fasta。
  2. 运行 `gtdbtk classify_wf`；输出 `gtdbtk.bac120.summary.tsv`。
* **工具**：`gtdbtk_classify(genome_dir, gtdb_ref="R220")`
* **Refdata**

  * 锁定 GTDB 参考库版本（Rxxx）；冻结 summary.tsv。
  * 评测：以 `assembly_id` 为键对比 `classification` 字段 exact match。

---

## D. 数据注释 / 分析类

### 9) InterProScan 注释（蓝细菌蛋白）

* **步骤**

  1. 输入 proteome FASTA。
  2. `interproscan.sh -i proteome.fa -f tsv`，得到 domain/GO/Pathway 注释。
* **工具**：`interproscan_annotate(proteome_fasta, formats=["tsv"])`
* **Refdata**

  * 冻结 `*.tsv`（排序规则：`protein_id, start, end, signature`）。
  * 评测：键列（`protein_id+signature+start+end`）exact match；允许数据库小版本差异时使用 `Jaccard≥0.98`。

---

### 10) KEGG 富集（COVID-19 DEGs）

* **步骤**

  1. 输入 DEG 列表（Gene Symbol + 上/下调）。
  2. 用 `clusterProfiler`（R）或 `gseapy`（Py）进行 pathway enrichment。
* **工具**：`kegg_enrich(deg_list, organism="hsa")`
* **Refdata**

  * 固定 KEGG 版本/日期与 p 值修正方法（BH）；冻结富集表（排序 by `p.adjust`）。
  * 评测：前 N 条路径名集合一致；富集统计值±1e-8 容差。

---

### 11) ARTS 检测耐药相关 BGC（RefSeq）

* **步骤**

  1. 输入细菌基因组。
  2. 运行 ARTS（或其 Docker）；输出 BGC 列表/坐标。
* **工具**：`arts_detect_bgc(genome_dir)`
* **Refdata**

  * 冻结 BGC 注释 TSV（规范列：`genome, region_id, product, start, end`）。
  * 评测：以 `genome+region_id` 为键 exact match；允许 product 字段小差异（可取集合相等）。

---

## E. 数据比对 / 关联类

### 12) 蛋白 → 基因组比对（Exonerate protein2genome）

* **步骤**

  1. 输入 RefSeq 蛋白 FASTA + Ensembl 基因组。
  2. `exonerate --model protein2genome prot.fa genome.fa --showvulgar no --showalignment no --ryo "<query>\t<target>\t%pi\t%qal\t%qae\t%tal\t%tae\n"`
* **工具**：`exonerate_protein2genome(protein_fa, genome_fa)`
* **Refdata**

  * 冻结对齐结果 TSV（列：`q_id, t_id, pct_identity, q_start, q_end, t_start, t_end`）。
  * 评测：以 `q_id+t_id+q_start+q_end` 为键 exact match；容许 `pct_identity`±0.01。

---

### 13) 基因–药物靶点映射（DrugCentral）

* **步骤**

  1. 输入 DEG 列表（Symbol）。
  2. 使用 DrugCentral 本地映射表（或 API）获得基因→药物列表；保存 `gene\tdrug\tmechanism`。
* **工具**：`drugcentral_map_targets(genes, local_tsv=...)`
* **Refdata**

  * 冻结一版 DrugCentral 映射表；冻结映射输出。
  * 评测：键列 `gene+drug` exact match；`mechanism` 作为补充字段（可忽略大小写）。

---

### 14) 宿主–性状匹配（AMDB × EltonTraits）

* **步骤**

  1. 输入宏基因组样本宿主（学名）与 EltonTraits（物种性状表）。
  2. 先用 NCBI Taxonomy 规范化学名（属/种对齐），再做 join。
* **工具**：`host_trait_join(hosts_tsv, traits_tsv, taxnorm_db=ncbi_taxdump)`
* **Refdata**

  * 冻结 join 后表（规范列：`sample_id, host, trait1, trait2, ...`）。
  * 评测：行键 `sample_id` 集合一致；缺失率（NA占比）≤ 参考值 + 0.5%。

---

## F. 数据检索 / 整理类

### 15) BioProject 关键词检索（“human gut microbiome”）

* **步骤**

  1. `esearch -db bioproject -query "human gut microbiome"` 获取 BioProject IDs。
  2. `esummary` 拉取标题/物种信息，保存 `bioproject_id\ttitle\ttaxon`。
* **工具**：`bioproject_keyword_search(query, date_from, date_to)`
* **Refdata**

  * 冻结 ID 列表与元信息（注意日期范围）；保存表。
  * 评测：ID 集合 Jaccard≥0.98（允许非常小的时间漂移）；标题/物种作为附加核对。

---

# 参考：工具封装接口（示例草图）

> 放在 `tools/*.py`，每个工具职责单一、输入输出路径明确。这里给一两个代表示例，其他同理扩展。

**示例：序列长度过滤**

```python
# tools/seq_length_filter.py
from Bio import SeqIO

def seq_length_filter(input_fasta, output_fasta, output_ids, min_len=50, max_len=500):
    kept = []
    with open(output_fasta, "w") as fo:
        for rec in SeqIO.parse(input_fasta, "fasta"):
            L = len(rec.seq)
            if min_len <= L <= max_len:
                kept.append(rec.id)
                SeqIO.write(rec, fo, "fasta")
    with open(output_ids, "w") as f:
        f.write("\n".join(sorted(kept)) + "\n")
    return output_fasta, output_ids
```

**示例：Exonerate 封装**

```python
# tools/exonerate_protein2genome.py
import subprocess, shlex

def exonerate_protein2genome(protein_fa, genome_fa, output_tsv):
    cmd = (
        'exonerate --model protein2genome '
        f'{shlex.quote(protein_fa)} {shlex.quote(genome_fa)} '
        '--showvulgar no --showalignment no '
        '--ryo "%qi\\t%ti\\t%pi\\t%qal\\t%qae\\t%tal\\t%tae\\n"'
    )
    with open(output_tsv, "w") as fo:
        subprocess.run(cmd, shell=True, check=True, stdout=fo)
    return output_tsv
```

> 评测脚本（`eval/*.py`）统一读取运行结果与 `benchmark/refdata/<task_id>/` 下的参考产物，按本任务的“评测”要求做比较并输出 PASS/FAIL 与差异统计。

---

# 如何准备参考答案（Refdata）的通用流程

1. **固定输入与版本**

   * 选定**小而代表性**的样本集；记录来源、下载日期、工具/库版本（写入 `meta/versions.json`）。
2. **一次性跑通**

   * 用“黄金环境”（常用 Docker/Conda 锁定）运行工具，产出结果。
3. **标准化与瘦身**

   * 对表格：排序（稳定键）、去重、统一列顺序与小数格式；对 FASTA：只保留目标子集；对对齐结果：只保留关键列。
   * 对下载类：从大文件抽出 **manifest** 清单（ID+MD5），减少存储。
4. **冻结**

   * 放入 `benchmark/refdata/<task_id>/`；保存 `README.md` 说明输入、命令、参数、版本与统计（行数、哈希）。
5. **定义评测规则**

   * 在 `eval/<task_id>.py` 固化该任务的比较逻辑（见各任务的“评测”条目）。
6. **回归测试**

   * 每次工具或依赖更新，先在 CI 运行**小样本**（smoke test），如通过再考虑更新 refdata（并升版本号）。

---

# 你可以直接落地的下一步

1. **我可以把这 15 个任务**整理成一份 `tasks.yaml`（每行一个任务，含：Category、Task\_type、工具名、输入/输出约定、评测器名），再给出 15 个 `eval_*.py` 骨架。
2. **我可以起草一个最小 agent**（Python）：

   * 读取 `tasks.yaml` → 调用 `tools/` 中的封装 → 把产物放到 `benchmark/results/<task_id>/run-<timestamp>/`
   * 然后自动调用对应 `eval/<task_id>.py` 与 `refdata/<task_id>/` 比较并给出报告。

如果你愿意，我就按上面的规范把 **`tasks.yaml` + 15 个工具函数骨架 + 评测脚本骨架** 生成出来，直接放到一个可运行的最小仓库结构里（你只需把大工具装好/映射好路径即可）。
