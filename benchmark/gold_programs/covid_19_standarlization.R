# 1. 安装必要的包
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")

if (!requireNamespace("Biobase", quietly = TRUE))
    BiocManager::install("Biobase", update = FALSE, ask = FALSE, repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")

if (!requireNamespace("pheatmap", quietly = TRUE))
    install.packages("pheatmap", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")

# 2. 加载包和数据
library(Biobase)
library(pheatmap)

# 加载 Bioconductor 自带的示例数据集 (sample.ExpressionSet)
data(sample.ExpressionSet)

# 3. 数据处理：提取表达矩阵
expr_matrix <- exprs(sample.ExpressionSet)

# 我们可以先看一眼数据长什么样
print(paste("基因数量:", nrow(expr_matrix)))
print(paste("样本数量:", ncol(expr_matrix)))

# 4. 绘图：选取前 30 个差异最大的基因画热图
# (计算方差，取方差最大的30个)
var_genes <- apply(expr_matrix, 1, var)
top30_genes <- names(sort(var_genes, decreasing = TRUE))[1:30]
mat_to_plot <- expr_matrix[top30_genes, ]

# 绘制热图
pheatmap(mat_to_plot, 
         scale = "row", # 按行归一化
         main = "Hello Bioinfo: Gene Expression Heatmap",
         show_rownames = TRUE)

print("Bioconductor 环境配置完美！")