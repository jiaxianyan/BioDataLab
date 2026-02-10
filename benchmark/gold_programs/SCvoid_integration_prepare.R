# BiocManager::install(c("GEOquery"))
# install.packages(
#   "GEOquery_2.70.0.tar.gz",
#   repos = NULL,
#   type = "source"
# )
library(GEOquery)
library(Matrix)
library(Seurat)
library(data.table)

gse_id <- "GSE182298"  # <- 换成你的
outdir <- file.path('benchmark/dataset/SCovid', gse_id)
outdir
# # dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# # 下载补充文件到 outdir
# # getGEOSuppFiles(GEO = gse_id, makeDirectory = FALSE, baseDir = outdir)

# # 找到可能的 10x 目录：含 matrix.mtx(.gz) 的目录
# mtx_files <- list.files(outdir, pattern = "matrix\\.mtx(\\.gz)?$", recursive = TRUE, full.names = TRUE)
# tenx_dirs <- unique(dirname(mtx_files))
# tenx_dirs

process_gse_flat_10x <- function(outdir,
                                 method = c("symlink", "copy"),
                                 overwrite = FALSE) {
  method <- match.arg(method)

  # 1) 找到所有 matrix 文件（用你这种命名规律）
  mtx <- list.files(outdir, pattern = "-matrix\\.mtx\\.gz$", full.names = TRUE)
  if (length(mtx) == 0) stop("没找到 *-matrix.mtx.gz 文件：", outdir)

  # 2) 推断前缀（样本ID）
  prefix <- sub("-matrix\\.mtx\\.gz$", "", basename(mtx))
  bc  <- file.path(outdir, paste0(prefix, "-barcodes.tsv.gz"))
  feat <- file.path(outdir, paste0(prefix, "-features.tsv.gz"))

  ok <- file.exists(bc) & file.exists(feat)
  if (!all(ok)) {
    message("有些样本三件套不齐，列出来：")
    print(data.frame(sample = prefix[!ok],
                     matrix = basename(mtx[!ok]),
                     barcodes_exists = file.exists(bc[!ok]),
                     features_exists = file.exists(feat[!ok])))
  }

  prefix <- prefix[ok]
  mtx <- mtx[ok]
  bc <- bc[ok]
  feat <- feat[ok]

  # 3) 为每个样本建目录并写入标准 10x 文件名
  created <- 0L
  for (i in seq_along(prefix)) {
    sample <- prefix[i]
    sample_dir <- file.path(outdir, sample)
    dir.create(sample_dir, showWarnings = FALSE, recursive = TRUE)

    target_mtx  <- file.path(sample_dir, "matrix.mtx.gz")
    target_bc   <- file.path(sample_dir, "barcodes.tsv.gz")
    target_feat <- file.path(sample_dir, "features.tsv.gz")

    targets <- c(target_mtx, target_bc, target_feat)
    if (!overwrite && any(file.exists(targets))) {
      message("跳过（已存在）：", sample_dir)
      next
    }

    # 如果 overwrite=TRUE，先删掉旧的目标文件（包括旧链接）
    if (overwrite) {
      for (t in targets) if (file.exists(t) || file.symlink(t)) unlink(t)
    }

    if (method == "symlink") {
      # 软链接：不占空间
      file.symlink(from = normalizePath(mtx[i]), to = target_mtx)
      file.symlink(from = normalizePath(bc[i]),  to = target_bc)
      file.symlink(from = normalizePath(feat[i]), to = target_feat)
    } else {
      # 拷贝：占空间
      file.copy(mtx[i], target_mtx, overwrite = TRUE)
      file.copy(bc[i],  target_bc,  overwrite = TRUE)
      file.copy(feat[i], target_feat, overwrite = TRUE)
    }

    created <- created + 1L
  }

  message("完成：成功处理样本数 = ", created,
          "（method=", method, ", overwrite=", overwrite, "）")
  invisible(TRUE)
}

# 用法：
outdir <- "benchmark/dataset/SCovid/GSE182298"
process_gse_flat_10x(outdir, method = "symlink", overwrite = FALSE)

