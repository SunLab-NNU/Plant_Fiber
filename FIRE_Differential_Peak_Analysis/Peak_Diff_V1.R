#!/usr/bin/env Rscript

suppressWarnings(suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(DESeq2)
  library(viridis)
  library(knitr)
  library(crayon)
  library(data.table)
  library(dendsort)
  library(dplyr)
  library(edgeR)
  library(extrafont)
  library(ggplot2)
  library(ggpubr)
  library(ggsci)
  library(gmodels)
  library(readr)
  library(rtracklayer)
  library(scatterplot3d)
  library(stringr)
  library(tidyverse)
  library(tximport)
  library(ggrepel)
  library(magrittr)
  library(clusterProfiler)
  library(biomaRt)
  library(org.At.tair.db)
  library(RColorBrewer)
  library(circlize)
  library(scales)
}))

suppressMessages(loadfonts())


setwd(getwd()) 
source("./RNA_Seq_Ara/theme_linhua.R")

theme_set(theme_linhua(base_size = 18, base_family = "", legend = "right"))
theme_update(
  panel.grid.minor.x = element_blank(),
  panel.grid.minor.y = element_blank(),
  panel.grid.major.x = element_blank(),
  panel.grid.major.y = element_blank(),
  panel.spacing.x = unit(0.5, "lines"),
  panel.spacing.y = unit(0.5, "lines"),
  strip.background = element_blank(),
  panel.border = element_rect(colour = "black"),
  aspect.ratio = 1
)


args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript script.R <input_file> <control_name> <treat_name>\n
       Example: Rscript script.R peaks.txt Control Treat")
}

input_file <- args[1]
A <- args[2] 
B <- args[3]


run_analysis <- function(input_file, A, B) {

  heat_fc <- fread(input_file, select = c(1, seq(2, 7)), skip = 1)
  colnames(heat_fc) <- c("seqnames", "start", "end", paste0(A, "-1"), paste0(A, "-2"), paste0(B, "-1"), paste0(B, "-2"))
  
  heat_fc$ID <- paste0("peak", 1:nrow(heat_fc))
  heat_fc_df <- as.data.frame(heat_fc)
  

  count_data <- heat_fc_df[, 4:7]
  rownames(count_data) <- heat_fc$ID
  
  count_data[] <- lapply(count_data, function(x) { 
    x[is.nan(x)] <- 0 
    return(x) 
  })
  

  sample_info <- data.frame(
    row.names = colnames(count_data),
    condition = factor(rep(c(A, B), each = 2), levels = c(A, B))
  )
  

  dds <- DESeqDataSetFromMatrix(
    countData = round(count_data), 
    colData = sample_info, 
    design = ~ condition
  )
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("condition", B, A))
  

  T_VS_C <- as.data.table(as.data.frame(res), keep.rownames = "ID")[
    , "direction" := "NoChange"
  ][
    , "Condition" := paste0(A, "_vs_", B)
  ]
  
  merged_data <- left_join(T_VS_C, heat_fc_df, by = "ID") %>% as.data.table()

for (C in seq(0.05, 1, by = 0.05)) {
  merged_data[, direction := NA_character_] 
  merged_data[pvalue < C & log2FoldChange > 1, direction := "Up"]
  merged_data[pvalue < C & log2FoldChange < -1, direction := "Down"]

  output_prefix <- paste0(A, "_VS_", B)

  up_peaks <- merged_data[direction == "Up", .(seqnames, start, end, ID, log2FoldChange)]
  down_peaks <- merged_data[direction == "Down", .(seqnames, start, end, ID, log2FoldChange)]

  C_label <- formatC(C, format = "f", digits = 2)

dir.create("Pvalue_C_DiffPeak", showWarnings = FALSE)
fwrite(up_peaks, file.path("Pvalue_C_DiffPeak", paste0("P", C_label, "_", output_prefix, "_up_Peaks.bed")), sep = "\t", col.names = FALSE)
fwrite(down_peaks, file.path("Pvalue_C_DiffPeak", paste0("P", C_label, "_", output_prefix, "_down_Peaks.bed")), sep = "\t", col.names = FALSE)

}

}


run_analysis(input_file, A, B)
