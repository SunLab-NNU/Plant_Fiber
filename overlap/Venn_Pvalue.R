library(GenomicRanges)
library(VennDiagram)
library(grid)
library(rtracklayer)
library(regioneR)
library(optparse) 
# 定义命令行参数解析
option_list <- list(
  make_option(c("-a", "--A_bed"), type = "character", help = "Path to first BED file (A)"),
  make_option(c("-b", "--B_bed"), type = "character", help = "Path to second BED file (B)"),
  make_option(c("-g", "--genome"), type = "character", help = "Path to genome length file"),
  make_option(c("-o", "--output"), type = "character", default = "Venn_with_Pvalue.pdf", 
              help = "Output PDF filename [default: %default]")
)

# 解析参数
opt <- parse_args(OptionParser(option_list = option_list))

# 检查必要参数
if (is.null(opt$A_bed) || is.null(opt$B_bed) || is.null(opt$genome)) {
  stop("Required arguments missing! Usage:
       Rscript script.R -a <A.bed> -b <B.bed> -g <genome.len> [-o output.pdf]", call. = FALSE)
}

# 读取 BED 文件
A <- import(opt$A_bed)
B <- import(opt$B_bed)

# 读取基因组长度文件并构建 GRanges
genome_data <- read.table(opt$genome, header = FALSE, stringsAsFactors = FALSE, quote="")
colnames(genome_data) <- c("seqnames", "length")

genome <- GRanges(
  seqnames = genome_data$seqnames,
  ranges = IRanges(start = 1, end = genome_data$length)
)

# 运行 permutation test
pt <- permTest(
  A = A,
  B = B,
  ntimes = 1000,
  randomize.function = randomizeRegions,
  evaluate.function = numOverlaps,
  genome = genome,
  force.parallel = FALSE
)

# 计算统计量
overlap_count <- length(findOverlaps(A, B, minoverlap = 1))
p_val <- signif(pt$numOverlaps$pval, 3)
total_A <- length(A)
total_B <- length(B)

# 绘制 Venn 图并保存
pdf(opt$output, width = 10, height = 10)
venn.plot <- draw.pairwise.venn(
  area1 = total_A,
  area2 = total_B,
  cross.area = overlap_count,
  category = c(basename(opt$A_bed), basename(opt$B_bed)),  # 使用文件名作为标签
  fill = c("skyblue", "salmon"),
  alpha = 0.5,
  lty = "solid",
  cex = c(2, 2, 2),
  label.col = c("black", "black", "red"),
  cat.cex = 1.5,
  cat.pos = c(-30, 30),
  cat.dist = 0.05,
  margin = 0.1
)
grid.text(
 label = paste0("Overlap = ", overlap_count, "\n", "p = ", formatC(p_val, format = "e", digits = 2)),
  x = 0.5, y = 0.8,
  gp = gpar(fontsize = 14, col = "black")
)
dev.off()

message(paste0("Venn diagram saved to: ", normalizePath(opt$output)))
