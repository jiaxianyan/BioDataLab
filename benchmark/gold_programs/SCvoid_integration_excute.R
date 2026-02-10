library(Seurat)
library(Matrix)

outdir <- "benchmark/dataset/SCovid/GSE182298"
samples <- list.dirs(outdir, recursive = FALSE, full.names = FALSE)

# 组织猜测：按你的样本名里真实包含的关键词来写
guess_tissue <- function(sample_name) {
  if (grepl("Olfactory_Bulb|嗅球", sample_name, ignore.case = TRUE)) return("Olfactory_Bulb")
  if (grepl("Frontal_Lobe|额叶", sample_name, ignore.case = TRUE))   return("Frontal_Lobe")
  if (grepl("liver|肝", sample_name, ignore.case = TRUE))            return("Liver")
  return("Other")
}

sample_tissue <- data.frame(
  sample = samples,
  tissue = vapply(samples, guess_tissue, character(1)),
  stringsAsFactors = FALSE
)

print(sample_tissue)

process_one_sample <- function(sample_dir, sample_name, tissue,
                               min_features = 200,
                               min_cells_per_gene = 3,
                               mt_liver = 50,
                               mt_other = 20) {

  mtx <- Read10X(sample_dir)
  obj <- CreateSeuratObject(mtx, project = sample_name)

  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

  mt_thresh <- ifelse(tissue %in% c("Liver", "liver", "Hepatocyte", "肝"),
                      mt_liver, mt_other)

  obj <- subset(obj, subset = nFeature_RNA >= min_features & percent.mt <= mt_thresh)

  counts <- GetAssayData(obj, layer = "counts")
  keep_genes <- Matrix::rowSums(counts > 0) >= min_cells_per_gene
  obj <- obj[keep_genes, ]

  obj$tissue <- tissue
  obj$mt_threshold <- mt_thresh
  obj
}

# ✅ 关键：在这里真正创建 objs
objs <- setNames(vector("list", nrow(sample_tissue)), sample_tissue$sample)

for (i in seq_len(nrow(sample_tissue))) {
  s <- sample_tissue$sample[i]
  t <- sample_tissue$tissue[i]
  sample_dir <- file.path(outdir, s)

  objs[[s]] <- process_one_sample(
    sample_dir = sample_dir,
    sample_name = s,
    tissue = t
  )
}

# ✅ 给每个样本的 cell 加前缀，避免 merge 冲突
objs <- lapply(names(objs), function(nm) RenameCells(objs[[nm]], add.cell.id = nm))
names(objs) <- sample_tissue$sample

# ✅ 合并
combined <- Reduce(function(a, b) merge(a, b), objs)

combined

passed <- colnames(combined)

# 提取 GSM 号
gsm <- sub("^(GSM[0-9]+).*", "\\1", passed)

# 提取 barcode（最后一个 "_" 后）
barcode <- sub(".*_", "", passed)

# 拼成 GSM_barcode
passed_gsm_barcode <- paste0(gsm, "_", barcode)

# 导出
write.table(
  passed_gsm_barcode,
  file = file.path(outdir, "qc_passed_barcodes.csv"),
  sep = ",",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)
