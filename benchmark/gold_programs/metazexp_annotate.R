setwd("~/biodatalab/benchmark/dataset/MetazExp/count")

files <- c(
  "DRX198149.gene.expr.tsv",
  "DRX198150.gene.expr.tsv",
  "DRX198151.gene.expr.tsv",
  "DRX198152.gene.expr.tsv",
  "DRX198153.gene.expr.tsv",
  "DRX198154.gene.expr.tsv"
)

# 读入并只保留 geneId + readCount
lst <- lapply(files, function(f){
  df <- read.delim(f, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  df <- df[, c("geneId","readCount")]
  colnames(df)[2] <- sub("\\.gene\\.expr\\.tsv$", "", f)  # 列名改成 SRX...
  df
})

# 按 geneId 合并（全连接），缺失补 0
merged <- Reduce(function(x,y) merge(x,y, by="geneId", all=TRUE), lst)
merged[is.na(merged)] <- 0

# geneId 作为行名，剩下都是整数 counts
count_mat <- as.matrix(merged[,-1])
rownames(count_mat) <- merged$geneId
storage.mode(count_mat) <- "integer"

# 保存给 DESeq2
write.table(count_mat, file="counts_matrix.tsv", sep="\t", quote=FALSE, col.names=NA)

coldata <- data.frame(
  sample = c(
    "DRX198149","DRX198150","DRX198151",
    "DRX198152","DRX198153","DRX198154"
  ),
  condition = c(
    "cancer","cancer","cancer",
    "control","control","control"
  )
)

rownames(coldata) <- coldata$sample


if (!requireNamespace("DESeq2", quietly = TRUE))
    BiocManager::install("DESeq2")
library(DESeq2)

dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = coldata,
  design = ~ condition
)

dds <- DESeq(dds)

res <- results(dds, contrast = c("condition","cancer","control"))
res <- res[order(res$padj), ]
head(res)

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[, c(
  "gene_id",
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "stat",
  "pvalue",
  "padj"
)]

write.table(
  res_df,
  file = "deseq2_results_all_genes.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
