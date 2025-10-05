#!/usr/bin/env Rscript

library(VennDiagram)
library(grid)

args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 2){
  stop("请提供两个 TXT 文件路径作为参数，例如：Rscript script.R file1.txt file2.txt")
}

file_A <- args[1]
file_B <- args[2]

# 读取两个 TXT 文件，假设每个文件只有一列基因 ID
A <- unique(read.table(file_A, stringsAsFactors = FALSE)[,1])
B <- unique(read.table(file_B, stringsAsFactors = FALSE)[,1])

# 设置拟南芥基因总量
genome_size <- 27000  # 可以根据实际基因数调整

# 计算交集和各集合大小
overlap_count <- length(intersect(A, B))
total_A <- length(A)
total_B <- length(B)

overlap_genes <- intersect(A, B)
print(overlap_genes)
# 超几何检验 P 值
p_val <- phyper(overlap_count - 1, total_A, genome_size - total_A, total_B, lower.tail = FALSE)
p_val <- signif(p_val, 3)

# 生成输出文件名
get_filename_without_ext <- function(filepath) {
  fname <- basename(filepath)
  sub("\\.[^.]*$", "", fname)
}

fileA_name <- get_filename_without_ext(file_A)
fileB_name <- get_filename_without_ext(file_B)
pdf_file <- paste0(fileA_name, "_vs_", fileB_name, "_Venn_with_Pvalue.pdf")

# 保存 PDF
pdf(pdf_file, width = 10, height = 10)

draw.pairwise.venn(
  area1 = total_A,
  area2 = total_B,
  cross.area = overlap_count,
  category = c(fileA_name, fileB_name),
  fill = c("skyblue", "salmon"),
  alpha = 0.5,
  lty = "solid",
  cex = c(2, 2, 2),
  label.col = c("black", "black", "black"),
  cat.cex = 1.5,
  cat.pos = c(-30, 30),
  cat.dist = 0.05,
  margin = 0.1
)

grid.text(
  label = paste0("overlap = ", overlap_count, "\n", "p = ", formatC(p_val, format = "e", digits = 2)),
  x = 0.5, y = 0.8,
  gp = gpar(fontsize = 14, col = "black")
)

dev.off()

